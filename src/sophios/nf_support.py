"""What the Nextflow backend does with each field it considers: one table.

Every row gives a field one of the three states of the design's closed-world
contract (design_docs/sophios_nextflow_backend.md §3) and a one-line note.

LOWERED: a lowering consumes it.
COMPILE_TIME: Sophios uses it while compiling; the Nextflow output carries nothing for it.
REJECTED: a finding at the field, whose message is the row's note, and no artifacts.

`FIELDS` is keyed by the `cwl_utils.parser.cwl_v1_2` class whose fields each part classifies.
Capability analysis (`sophios.utils_nf`) consumes exactly the fields that are not rejected, and
docs/nextflow.md renders this table (docs/conf.py). `tests/core/nextflow/test_support_table.py`
holds `FIELDS` to the fields cwl_utils declares and to Sophios's own `SUPPORT_MATRIX`.
"""
from dataclasses import dataclass
from enum import StrEnum
from types import MappingProxyType
from typing import Final, Mapping


class State(StrEnum):
    """What the Nextflow backend does with one field."""

    LOWERED = 'lowered'
    COMPILE_TIME = 'compile-time'
    REJECTED = 'rejected'


@dataclass(frozen=True, slots=True)
class Row:
    """One field's state and what it means for that field."""

    state: State
    note: str


_L, _C, _R = State.LOWERED, State.COMPILE_TIME, State.REJECTED

_DOC = Row(_C, 'documentation')
_INTENT = Row(_C, 'documentation: an ontology term for what it is for')
_FORMAT = Row(_C, 'read by EDAM inference')
_STREAMABLE = Row(_C, 'a hint to the runner; a task reads a staged file')
_SECONDARY_FILES = Row(_R, 'secondaryFiles are not supported: stage each file as its own input')
_LOAD_CONTENTS = Row(_R, "loadContents on an input is not supported: an input's contents would have to "
                         'become a value channel read at staging time')
_LOAD_LISTING = Row(_R, 'loadListing is not supported: a Directory is staged whole')

FIELDS: Final[Mapping[str, Mapping[str, Row]]] = MappingProxyType({
    'Workflow': MappingProxyType({
        'class': Row(_L, '`Workflow`'),
        'cwlVersion': Row(_L, '`v1.2`, whose semantics the lowering follows'),
        'doc': _DOC,
        'hints': Row(_R, 'workflow hints are not supported: put each one on the tool that needs it'),
        'id': Row(_L, 'names the Nextflow workflow'),
        'inputs': Row(_L, 'workflow parameters'),
        'intent': _INTENT,
        'label': _DOC,
        'outputs': Row(_L, 'workflow emits'),
        'requirements': Row(_L, 'feature declarations, each inert; see the workflow requirement classes'),
        'steps': Row(_L, "one process per tool step; a nested workflow's steps are inlined"),
    }),
    'WorkflowStep': MappingProxyType({
        'doc': _DOC,
        'hints': Row(_R, "a step's hints are not supported: put them on the step's tool"),
        'id': Row(_L, 'names the process'),
        'in': Row(_L, "the process's input channels"),
        'label': _DOC,
        'out': Row(_L, "the process's emits"),
        'requirements': Row(_R, "a step's requirements are not supported: put them on the step's tool"),
        'run': Row(_L, 'a `CommandLineTool`, or a nested workflow, inlined'),
        'scatter': Row(_L, 'one task per invocation, indexed so a gather keeps order'),
        'scatterMethod': Row(_L, '`dotproduct`, `flat_crossproduct` or `nested_crossproduct`'),
        'when': Row(_L, 'the safe JavaScript subset, evaluated per invocation'),
    }),
    'WorkflowStepInput': MappingProxyType({
        'default': Row(_R, "a step input's default is not supported"),
        'id': Row(_L, 'the tool input it binds'),
        'label': _DOC,
        'linkMerge': Row(_R, "a step input's linkMerge is not supported"),
        'loadContents': Row(_R, "a step input's loadContents is not supported"),
        'loadListing': Row(_R, "a step input's loadListing is not supported"),
        'pickValue': Row(_R, "a step input's pickValue is not supported"),
        'source': Row(_L, 'one workflow input or step output'),
        'valueFrom': Row(_R, "a step input's valueFrom is not supported"),
    }),
    'WorkflowInputParameter': MappingProxyType({
        'default': Row(_L, "the parameter's value when none is given; rejected on a nested workflow's input"),
        'doc': _DOC,
        'format': _FORMAT,
        'id': Row(_L, 'names the parameter'),
        'inputBinding': Row(_R, "a workflow input's inputBinding is not supported"),
        'label': _DOC,
        'loadContents': Row(_R, "a workflow input's loadContents is not supported"),
        'loadListing': _LOAD_LISTING,
        'secondaryFiles': _SECONDARY_FILES,
        'streamable': _STREAMABLE,
        'type': Row(_L, "the parameter's channel: `path` for `File` and `Directory`, `val` otherwise"),
    }),
    'WorkflowOutputParameter': MappingProxyType({
        'doc': _DOC,
        'format': _FORMAT,
        'id': Row(_L, 'names the emit'),
        'label': _DOC,
        'linkMerge': Row(_R, "a workflow output's linkMerge is not supported"),
        'outputSource': Row(_L, 'one step output'),
        'pickValue': Row(_R, "a workflow output's pickValue is not supported"),
        'secondaryFiles': _SECONDARY_FILES,
        'streamable': _STREAMABLE,
        'type': Row(_L, "the emitted value's type"),
    }),
    'CommandLineTool': MappingProxyType({
        'arguments': Row(_L, 'command tokens, ordered as CWL orders them'),
        'baseCommand': Row(_L, 'the first command tokens, as literal text'),
        'class': Row(_L, '`CommandLineTool`; any other class is rejected'),
        'cwlVersion': Row(_L, '`v1.2`, whose semantics the lowering follows'),
        'doc': _DOC,
        'hints': Row(_L, 'read as requirements; see the tool requirement classes'),
        'id': Row(_L, 'names the tool in diagnostics'),
        'inputs': Row(_L, "the process's inputs"),
        'intent': _INTENT,
        'label': _DOC,
        'outputs': Row(_L, "the process's outputs"),
        'permanentFailCodes': Row(_R, 'permanentFailCodes is not supported: a nonzero exit fails the task'),
        'requirements': Row(_L, 'see the tool requirement classes'),
        'stderr': Row(_L, "the script's stderr redirection"),
        'stdin': Row(_L, "the script's stdin redirection"),
        'stdout': Row(_L, "the script's stdout redirection"),
        'successCodes': Row(_R, 'successCodes is not supported: a nonzero exit fails the task'),
        'temporaryFailCodes': Row(_R, 'temporaryFailCodes is not supported: a nonzero exit fails the task'),
    }),
    'CommandInputParameter': MappingProxyType({
        'default': Row(_L, 'a workflow parameter for an input no step input binds'),
        'doc': _DOC,
        'format': _FORMAT,
        'id': Row(_L, 'names the port'),
        'inputBinding': Row(_L, 'command tokens'),
        'label': _DOC,
        'loadContents': _LOAD_CONTENTS,
        'loadListing': _LOAD_LISTING,
        'secondaryFiles': _SECONDARY_FILES,
        'streamable': _STREAMABLE,
        'type': Row(_L, "the port's qualifier: `path` for `File` and `Directory`, "
                    '`val` for `int`, `float`, `string` and `boolean`'),
    }),
    'CommandLineBinding': MappingProxyType({
        'itemSeparator': Row(_R, 'itemSeparator is not supported: an array binding gives each element its own word'),
        'loadContents': _LOAD_CONTENTS,
        'position': Row(_L, 'orders the tokens'),
        'prefix': Row(_L, "literal text before the value, or a `boolean`'s flag"),
        'separate': Row(_L, '`false` joins the prefix and the value into one word'),
        'shellQuote': Row(_L, '`false` only on literal text under `ShellCommandRequirement`'),
        'valueFrom': Row(_L, 'literal text or `$(inputs.<name>)`, and on a `boolean` only its bare '
                             'self-reference; a computed `valueFrom` is rejected'),
    }),
    'CommandOutputParameter': MappingProxyType({
        'doc': _DOC,
        'format': _FORMAT,
        'id': Row(_L, 'names the emit'),
        'label': _DOC,
        'outputBinding': Row(_L, "the output's capture"),
        'secondaryFiles': _SECONDARY_FILES,
        'streamable': _STREAMABLE,
        'type': Row(_L, '`File` or `Directory` by glob, or `string` by file-text capture'),
    }),
    'CommandOutputBinding': MappingProxyType({
        'glob': Row(_L, 'a path glob of literal text and input references'),
        'loadContents': Row(_L, 'file-text capture, beside `outputEval: $(self[0].contents)`'),
        'loadListing': _LOAD_LISTING,
        'outputEval': Row(_L, '`$(self[0])`, or `$(self[0].contents)` beside `loadContents`, over a literal glob'),
    }),
    'DockerRequirement': MappingProxyType({
        'class': Row(_L, 'the `container` directive'),
        'dockerFile': Row(_R, 'dockerFile is not supported: use dockerPull'),
        'dockerImageId': Row(_L, 'the image'),
        'dockerImport': Row(_R, 'dockerImport is not supported: use dockerPull'),
        'dockerLoad': Row(_R, 'dockerLoad is not supported: use dockerPull'),
        'dockerOutputDirectory': Row(_R, 'dockerOutputDirectory is not supported: a task writes to its work directory'),
        'dockerPull': Row(_L, 'the image'),
    }),
    'InitialWorkDirRequirement': MappingProxyType({
        'class': Row(_L, "stages the inputs a listing names"),
        'listing': Row(_L, 'entries naming one `File` or `Directory` input, under its basename or a literal name'),
    }),
    'InlineJavascriptRequirement': MappingProxyType({
        'class': Row(_L, 'a declaration only; each expression is judged where it appears'),
        'expressionLib': Row(_R, 'expressionLib is not supported: the safe JavaScript subset has no functions'),
    }),
    'NetworkAccess': MappingProxyType({
        'class': Row(_C, 'a Nextflow task has network access'),
        'networkAccess': Row(_C, 'the literal `true` only; any other value is rejected'),
    }),
    'ResourceRequirement': MappingProxyType({
        'class': Row(_L, 'the `cpus` and `memory` directives'),
        'coresMax': Row(_L, '`cpus`, when `coresMin` is absent; a whole number'),
        'coresMin': Row(_L, '`cpus`; a whole number'),
        'outdirMax': Row(_R, 'outdirMax is not supported'),
        'outdirMin': Row(_R, 'outdirMin is not supported'),
        'ramMax': Row(_L, '`memory`, when `ramMin` is absent'),
        'ramMin': Row(_L, '`memory`'),
        'tmpdirMax': Row(_R, 'tmpdirMax is not supported'),
        'tmpdirMin': Row(_R, 'tmpdirMin is not supported'),
    }),
    'ScatterFeatureRequirement': MappingProxyType({
        'class': Row(_L, 'a declaration only; each scatter is lowered at its step'),
    }),
    'ShellCommandRequirement': MappingProxyType({
        'class': Row(_L, 'admits `shellQuote: false` on literal text'),
    }),
    'SubworkflowFeatureRequirement': MappingProxyType({
        'class': Row(_L, 'a declaration only; each nested workflow is inlined'),
    }),
})

_WORKFLOW_FEATURE = 'is a workflow requirement: declare it on the workflow'

#: The requirement and hint classes a CommandLineTool may name, and the
#: requirement classes a Workflow may name; any other class is rejected. A
#: class that is not rejected has its fields in `FIELDS`.
REQUIREMENT_CLASSES: Final[Mapping[str, Mapping[str, Row]]] = MappingProxyType({
    'CommandLineTool': MappingProxyType({
        'DockerRequirement': Row(_L, 'the `container` directive'),
        'InitialWorkDirRequirement': Row(_L, 'input staging'),
        'InlineJavascriptRequirement': Row(_L, 'a declaration only'),
        'LoadListingRequirement': Row(_R, 'LoadListingRequirement is not supported: a Directory is staged whole'),
        'NetworkAccess': Row(_C, '`networkAccess: true` only'),
        'ResourceRequirement': Row(_L, 'the `cpus` and `memory` directives'),
        'ScatterFeatureRequirement': Row(_R, f'ScatterFeatureRequirement {_WORKFLOW_FEATURE}'),
        'ShellCommandRequirement': Row(_L, 'admits `shellQuote: false` on literal text'),
        'SoftwareRequirement': Row(_R, 'SoftwareRequirement is not supported: use DockerRequirement'),
        'SubworkflowFeatureRequirement': Row(_R, f'SubworkflowFeatureRequirement {_WORKFLOW_FEATURE}'),
    }),
    'Workflow': MappingProxyType({
        'InlineJavascriptRequirement': Row(_L, 'a declaration only; Sophios adds it for every `when`'),
        'ScatterFeatureRequirement': Row(_L, 'a declaration only'),
        'SubworkflowFeatureRequirement': Row(_L, 'a declaration only'),
    }),
})


def consumed(rows: Mapping[str, Row]) -> frozenset[str]:
    """The names a row lowers or uses at compile time: every one but the rejected."""
    return frozenset(name for name, row in rows.items() if row.state is not State.REJECTED)
