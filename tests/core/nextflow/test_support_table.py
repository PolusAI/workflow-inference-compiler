"""The backend's support table is complete, and capability analysis reads it (design §3)."""

# pylint: disable=missing-function-docstring

from collections.abc import Callable
import copy
import re
from typing import Any

import pytest
from cwl_utils.parser import cwl_v1_2

from sophios.lang import SUPPORT_MATRIX
from sophios.nf_support import FIELDS, REQUIREMENT_CLASSES, State
from sophios.utils_nf import compiled_source_to_nextflow
from sophios.wic_types import Yaml

from .testkit import step, synthetic_source, tool, workflow_doc


@pytest.mark.fast
@pytest.mark.parametrize("cwl_class", sorted(SUPPORT_MATRIX))
def test_the_sophios_classes_have_the_support_matrix_fields(cwl_class: str) -> None:
    assert set(FIELDS[cwl_class]) == set(SUPPORT_MATRIX[cwl_class])


@pytest.mark.fast
@pytest.mark.parametrize("cwl_class", sorted(FIELDS))
def test_every_field_cwl_utils_declares_is_classified(cwl_class: str) -> None:
    assert set(FIELDS[cwl_class]) == set(getattr(cwl_v1_2, cwl_class).attrs)


@pytest.mark.fast
def test_every_requirement_class_is_a_cwl_class_and_an_admitted_one_has_its_fields() -> None:
    for classes in REQUIREMENT_CLASSES.values():
        for name, row in classes.items():
            assert isinstance(getattr(cwl_v1_2, name), type), name
            assert row.state is State.REJECTED or name in FIELDS, name


# Design §3 names the compile-time fields: these wherever CWL declares them, and NetworkAccess.
COMPILE_TIME = {"doc", "format", "intent", "label", "streamable"}


@pytest.mark.fast
def test_the_compile_time_fields_are_the_ones_the_design_names() -> None:
    compile_time = {
        (cwl_class, name) for cwl_class, rows in FIELDS.items() for name, row in rows.items()
        if row.state is State.COMPILE_TIME
    }
    assert compile_time == {
        (cwl_class, name) for cwl_class, rows in FIELDS.items() for name in rows
        if name in COMPILE_TIME or cwl_class == "NetworkAccess"
    }


def _documents() -> tuple[Yaml, Yaml]:
    run = tool(
        "T",
        inputs={"x": {"type": "string", "inputBinding": {"position": 1}}},
        outputs={"o": {"type": "File", "outputBinding": {"glob": "o.txt"}}},
        requirements={},
    )
    workflow = workflow_doc(
        [step("T", **{"in": {"x": {"source": "x"}}, "out": ["o"]})],
        inputs={"x": {"type": "string"}},
        outputs={"o": {"type": "File", "outputSource": "T/o"}},
    )
    return workflow, run


# Where each class's value sits in `_documents()`, and the path a finding names it by.
POSITIONS: dict[str, tuple[str, Callable[[Yaml, Yaml], Yaml]]] = {
    "Workflow": ("workflow", lambda workflow, _run: workflow),
    "WorkflowStep": ("steps[0]", lambda workflow, _run: workflow["steps"][0]),
    "WorkflowStepInput": ("steps[0].in.x", lambda workflow, _run: workflow["steps"][0]["in"]["x"]),
    "WorkflowInputParameter": ("workflow.inputs.x", lambda workflow, _run: workflow["inputs"]["x"]),
    "WorkflowOutputParameter": ("workflow.outputs.o", lambda workflow, _run: workflow["outputs"]["o"]),
    "CommandLineTool": ("steps[0].run", lambda _workflow, run: run),
    "CommandInputParameter": ("steps[0].run.inputs.x", lambda _workflow, run: run["inputs"]["x"]),
    "CommandLineBinding": (
        "steps[0].run.inputs.x.inputBinding", lambda _workflow, run: run["inputs"]["x"]["inputBinding"]
    ),
    "CommandOutputParameter": ("steps[0].run.outputs.o", lambda _workflow, run: run["outputs"]["o"]),
    "CommandOutputBinding": (
        "steps[0].run.outputs.o.outputBinding", lambda _workflow, run: run["outputs"]["o"]["outputBinding"]
    ),
}


def _position(cwl_class: str) -> tuple[str, Callable[[Yaml, Yaml], Yaml]]:
    """A row's position; a requirement class's fields sit in the tool's requirements."""
    if cwl_class in POSITIONS:
        return POSITIONS[cwl_class]
    return (
        f"steps[0].run.requirements.{cwl_class}",
        lambda _workflow, run: run["requirements"].setdefault(cwl_class, {}),
    )


def _lower(workflow: Yaml, run: Yaml, **values: Any) -> Any:
    return compiled_source_to_nextflow(synthetic_source(workflow, [run], workflow_inputs={"x": "hi", **values}))


REJECTED = [
    (cwl_class, name) for cwl_class, rows in FIELDS.items() for name, row in rows.items()
    if row.state is State.REJECTED
]


@pytest.mark.fast
@pytest.mark.parametrize(("cwl_class", "name"), REJECTED, ids=[f"{c}.{n}" for c, n in REJECTED])
def test_a_rejected_field_is_a_finding_whose_message_is_its_note(cwl_class: str, name: str) -> None:
    workflow, run = _documents()
    path, place = _position(cwl_class)
    place(workflow, run)[name] = "x"
    finding = re.escape(f"- {path}.{name}: {FIELDS[cwl_class][name].note}")
    with pytest.raises(ValueError, match=rf"{finding}(\n|\Z)"):
        _lower(workflow, run)


SAMPLES = {
    "doc": "documentation", "label": "a label", "intent": ["http://edamontology.org/operation_0004"],
    "format": "http://edamontology.org/format_1929", "streamable": True,
}


@pytest.mark.fast
def test_compile_time_fields_leave_the_executable_model_unchanged() -> None:
    workflow, run = _documents()
    # A File input, so its job value carries `format` too: a list, as a `.wic` file's `!ii` writes it.
    run["inputs"]["f"] = {"type": "File", "inputBinding": {"position": 3}}
    workflow["inputs"]["f"] = {"type": "File"}
    workflow["steps"][0]["in"]["f"] = {"source": "f"}
    plain = _lower(copy.deepcopy(workflow), copy.deepcopy(run), f={"class": "File", "path": "f.txt"})
    for cwl_class, (_path, place) in POSITIONS.items():
        place(workflow, run).update({name: SAMPLES[name] for name in COMPILE_TIME & set(FIELDS[cwl_class])})
    run["requirements"]["NetworkAccess"] = {"networkAccess": True}
    run["hints"] = [{"class": "NetworkAccess", "networkAccess": True}]
    value = {"class": "File", "location": "f.txt", "format": [SAMPLES["format"]]}
    assert _lower(workflow, run, f=value) == plain


@pytest.mark.fast
@pytest.mark.parametrize("section", ["requirements", "hints"])
@pytest.mark.parametrize("value", [False, "$(inputs.x)", None])
def test_network_access_other_than_the_literal_true_is_rejected_by_name(section: str, value: Any) -> None:
    workflow, run = _documents()
    run[section] = {"NetworkAccess": {} if value is None else {"networkAccess": value}}
    with pytest.raises(ValueError, match=re.escape(
        f"steps[0].run.{section}.NetworkAccess.networkAccess: only the literal true is supported"
    )):
        _lower(workflow, run)


@pytest.mark.fast
def test_a_rejected_requirement_class_is_a_finding_whose_message_is_its_note() -> None:
    workflow, run = _documents()
    run["requirements"] = {"SoftwareRequirement": {"packages": []}}
    note = REQUIREMENT_CLASSES["CommandLineTool"]["SoftwareRequirement"].note
    with pytest.raises(ValueError, match=re.escape(f"steps[0].run.requirements.SoftwareRequirement: {note}")):
        _lower(workflow, run)
