"""Channels and pipeline safety under installed, pinned Nextflow.

design_docs/sophios_nextflow_backend.md §5.2 and §6, Inputs and channels and
Topology; §7 for the stub every process carries.
"""

# pylint: disable=missing-function-docstring

from dataclasses import replace
from pathlib import Path
import subprocess
from typing import Any

import pytest

from sophios.api.python._workflow_runtime import compile_workflow_result
from sophios.api.python.tool_builder import CommandLineTool, Input, Inputs, Output, Outputs, cwl
from sophios.api.python.workflow import Step, Workflow
from sophios.input_output_nf import NF_STUB_MESSAGE, write_nextflow_artifacts
from sophios.utils_nf import compilation_result_source, compiled_source_to_nextflow

from .testkit import execute_nextflow, runtime_workflow


def _run(
    flow: Workflow, directory: Path, *arguments: str, params: dict[str, Any] | None = None
) -> tuple["subprocess.CompletedProcess[str]", dict[str, list[str]]]:
    """Lower ``flow`` with ``params`` as its workflow input values, run it, and read each workflow output.

    Each output is the text of every value it emits, a file's contents for a path.
    """
    source = compilation_result_source(compile_workflow_result(flow))
    workflow = compiled_source_to_nextflow(source if params is None else replace(source, params=params))
    write_nextflow_artifacts(workflow, directory)
    script = directory / "workflow.nf"
    emits = [connection.to_port for connection in workflow.connections if hasattr(connection, "to_port")
             and not hasattr(connection, "to_process")]
    sinks = "".join(
        f"\n    {workflow.name}.out.{emit}.toList().subscribe {{ list -> new File('{emit}.sink').text = "
        f"list.collect {{ it instanceof java.nio.file.Path ? it.text : it.toString() }}.join('\\u0001') }}\n"
        for emit in emits
    )
    script.write_text(script.read_text(encoding="utf-8")[:-2] + sinks + "}\n", encoding="utf-8")
    result = execute_nextflow(directory, *arguments)
    observed = {
        emit: (directory / f"{emit}.sink").read_text(encoding="utf-8").split("\x01")
        for emit in emits
        if (directory / f"{emit}.sink").exists()
    }
    return result, observed


@pytest.mark.nextflow
@pytest.mark.serial
def test_a_dry_run_stops_at_its_first_task_and_writes_nothing(tmp_path: Path) -> None:
    write_nextflow_artifacts(runtime_workflow(), tmp_path)
    result = execute_nextflow(tmp_path, "-stub-run")
    assert result.returncode == 1, result.stdout + result.stderr
    assert NF_STUB_MESSAGE in result.stdout + result.stderr
    assert not list((tmp_path / "work").rglob("message.txt"))
    assert not list((tmp_path / "work").rglob("copy.txt"))


def _scattered_and_whole(item: Any, read: str) -> Workflow:
    """One array input scattered at one step and taken whole at another: ``read`` prints a value."""
    each_tool = CommandLineTool(
        "each_value", Inputs(value=Input(item, position=1)), Outputs(line=Output(cwl.file, glob="each.txt"))
    ).base_command(read).stdout("each.txt")
    whole_tool = CommandLineTool(
        "all_values", Inputs(every=Input(cwl.array(item), position=1)), Outputs(line=Output(cwl.file, glob="all.txt"))
    ).base_command(read).stdout("all.txt")
    each, whole = Step(each_tool, step_name="each_value"), Step(whole_tool, step_name="all_values")
    flow = Workflow([each, whole], "shared")
    each.inputs.value = flow.inputs.shared.as_type({"type": "array", "items": item})
    whole.inputs.every = flow.inputs.shared
    each.scatter_on(each.inputs.value)
    flow.outputs.each, flow.outputs.whole = each.outputs.line, whole.outputs.line
    return flow


@pytest.mark.nextflow
@pytest.mark.serial
def test_one_string_array_is_scattered_at_one_step_and_taken_whole_at_another(tmp_path: Path) -> None:
    flow = _scattered_and_whole(cwl.string, "echo")
    result, observed = _run(flow, tmp_path, params={"shared": ["x", "y z"]})
    assert result.returncode == 0, result.stdout + result.stderr
    assert observed["each"] == ["x\n", "y z\n"]
    assert observed["whole"] == ["x y z\n"]


@pytest.mark.nextflow
@pytest.mark.serial
def test_one_file_array_is_scattered_at_one_step_and_taken_whole_at_another(tmp_path: Path) -> None:
    for name in ("a.txt", "b.txt"):
        (tmp_path / name).write_text(f"{name}\n", encoding="utf-8")
    flow = _scattered_and_whole(cwl.file, "cat")
    files = [str(tmp_path / "a.txt"), {"class": "File", "path": str(tmp_path / "b.txt")}]
    result, observed = _run(flow, tmp_path, params={"shared": files})
    assert result.returncode == 0, result.stdout + result.stderr
    assert observed["each"] == ["a.txt\n", "b.txt\n"]
    assert observed["whole"] == ["a.txt\nb.txt\n"]
