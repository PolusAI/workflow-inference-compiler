"""Compiled CWL is plain CWL: no `sophios` key or prefix, and EDAM only where a format says `edam:`."""
from typing import Any

import pytest

from sophios import post_compile
from sophios.ir.artifacts import CompilationArtifact
from sophios.ir.emit import EDAM_NAMESPACE

from .hermetic import compile_hermetic, subworkflow_step


def _documents(node: CompilationArtifact) -> list[Any]:
    return [node.cwl, *(document for child in node.children for document in _documents(child))]


#: The fields a CWL document declares its prefixes and ontologies in.
DECLARATIONS = ('$namespaces', '$schemas')


def declared_beneath(artifact: CompilationArtifact) -> dict[str, Any]:
    """The `$namespaces` and `$schemas` of `artifact` and of every process beneath it, merged as
    embedding merges them: an embedded process may not declare its own, so they move up."""
    namespaces = dict(artifact.cwl.get('$namespaces', {}))
    schemas = list(artifact.cwl.get('$schemas', []))
    for child in artifact.children:
        below = declared_beneath(child)
        namespaces |= below.get('$namespaces', {})
        schemas += [schema for schema in below.get('$schemas', []) if schema not in schemas]
    return {field: value for field, value in zip(DECLARATIONS, (namespaces, schemas)) if value}


def _keys(value: Any) -> list[str]:
    """Every mapping key in `value`, at any depth."""
    match value:
        case dict():
            return [str(key) for key in value] + [key for item in value.values() for key in _keys(item)]
        case list():
            return [key for item in value for key in _keys(item)]
    return []


@pytest.mark.fast
@pytest.mark.parametrize('form', ['nested', 'embedded', 'flat', 'flat_embedded'])
def test_no_form_of_the_compiled_cwl_carries_a_sophios_key_or_prefix(form: str) -> None:
    """Whichever way the tree is written, nothing in it is Sophios's own."""
    nested = compile_hermetic({'steps': [subworkflow_step('child.wic', {'steps': [{'id': 'mk_file'}]}),
                                         {'id': 'count'}]}).artifact
    flat = post_compile.flatten_subworkflows(nested)
    artifact = {'nested': nested, 'embedded': post_compile.inline_artifact_runs(nested),
                'flat': flat, 'flat_embedded': post_compile.inline_artifact_runs(flat)}[form]
    keys = [key for document in _documents(artifact) for key in _keys(document)]
    assert not [key for key in keys if key.startswith('sophios')]
    assert not [document for document in _documents(artifact) if 'sophios' in document.get('$namespaces', {})]


@pytest.mark.fast
def test_a_workflow_no_format_of_says_edam_binds_no_prefix() -> None:
    """`count` takes an unformatted File and makes an int: nothing in the document is an EDAM CURIE, so
    there is no `$namespaces` or `$schemas` at all. (The corpus tests validate such documents with cwltool.)"""
    artifact = compile_hermetic({'steps': [{'id': 'count'}]}, 'plain').artifact
    assert '$namespaces' not in artifact.cwl and '$schemas' not in artifact.cwl


@pytest.mark.fast
@pytest.mark.parametrize('workflow', [
    {'inputs': {'f': {'type': 'File', 'format': 'edam:format_2330'}}, 'steps': [{'id': 'count', 'in': {'file': 'f'}}]},
    {'steps': [{'id': 'mk_file'}]},
    {'steps': [{'id': 'count', 'in': {'file': {'wic_inline_input': {
        'class': 'File', 'location': 'x.txt', 'format': 'edam:format_2330'}}}}]},
    {'inputs': {'f': {'type': 'File', 'default': {'class': 'File', 'location': 'x.txt', 'format': 'edam:format_2330'}}},
     'steps': [{'id': 'count', 'in': {'file': 'f'}}]},
    {'steps': [{'id': 'count', 'in': {'file': {'wic_raw_cwl': {
        'default': {'class': 'File', 'location': 'x.txt', 'format': 'edam:format_2330'}}}}}]},
], ids=['authored-input-format', 'promoted-output-format', 'job-file-format', 'input-default-format',
        'step-input-default-format'])
def test_a_workflow_a_format_of_says_edam_binds_the_prefix(workflow: dict[str, Any]) -> None:
    """An authored format, a format promoted from a tool, and a File in the job file: a runner reads each
    CURIE against the workflow's `$namespaces`, so each needs `edam` bound. No schema is added: an ontology
    is a dependency only its author adds."""
    cwl = compile_hermetic(workflow).artifact.cwl
    assert cwl['$namespaces'] == {EDAM_NAMESPACE[0]: EDAM_NAMESPACE[1]} and '$schemas' not in cwl
