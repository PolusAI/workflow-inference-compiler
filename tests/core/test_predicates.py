"""Claims about a single compilation.

These are the properties that are *not* equivalences. An equivalence
compares two compilations and needs a relation to say what may differ between
them; a predicate examines one compilation and needs none. Grouping them by
that fact rather than by which register row they came from is the point — a
reader can tell at a glance which properties compare two runs and which read
one.

One claim lives here, **made-name injectivity**: distinct authored steps and
ports never get the same made name. Edge soundness is `cwltool --validate`'s,
over the same generator (`test_emit.test_emit_validates_as_cwl_v1_2`); the
fixed-point guard's is `test_infer_phase.test_iteration_exhaustion_is_exactly_wic022`'s.
"""
import copy
import pytest
from hypothesis import example, given
from hypothesis import strategies as st

from sophios.ir.artifacts import CompilationArtifact
from sophios.utils_cwl import desugar_into_canonical_normal_form
from sophios.wic_types import Cwl, StepId, Tool, Tools, Yaml

from .hermetic import ORACLE, compile_hermetic, subworkflow_step
from .synthetic_tools import SYNTHETIC_NS, clt

# --------------------------------------------------------------------------
# Made-name injectivity
# --------------------------------------------------------------------------

#: Authored names (wic045 admits `[A-Za-z0-9_.-]`) over a small alphabet, so
#: that one name is often another's prefix.
_NAMES = st.from_regex(r'[ab][ab_.-]{0,3}', fullmatch=True)
_SEPARATORS = st.sampled_from('_.-')

#: A tool's inputs and outputs, each a set of authored port names.
_Signature = tuple[frozenset[str], frozenset[str]]


def _tools(signatures: dict[str, _Signature]) -> Tools:
    """A registry with one stub tool per stem, every port a string."""
    def tool(inputs: frozenset[str], outputs: frozenset[str]) -> Cwl:
        return clt({port: {'type': 'string', 'inputBinding': {'position': 1}} for port in inputs},
                   {port: {'type': 'string', 'outputBinding': {'outputEval': '$("x")'}} for port in outputs},
                   javascript=True)
    return {StepId(stem, SYNTHETIC_NS): Tool(f'/synthetic/{stem}.cwl',
                                             desugar_into_canonical_normal_form(tool(*signature)))
            for stem, signature in signatures.items()}


def _workflow(stems: list[str], signatures: dict[str, _Signature], calls: list[str],
              bodies: dict[str, Yaml]) -> Yaml:
    """Tool steps, each input bound to a literal, then a call of each subworkflow named in `calls`."""
    steps: list[Yaml] = [{'id': stem, 'in': {port: {'wic_inline_input': f'{stem}/{port}'}
                                             for port in sorted(signatures[stem][0])}}
                         for stem in stems]
    return {'steps': steps + [subworkflow_step(f'{name}.wic', copy.deepcopy(bodies[name])) for name in calls]}


@st.composite
def _signatures(draw: st.DrawFn) -> dict[str, _Signature]:
    """Tools whose names are built to meet: step `x` with port `y·z` beside
    step `x·y` with port `z`, `·` being `_`, `.` or `-`."""
    ports = st.frozensets(_NAMES, max_size=2)
    drawn = draw(st.dictionaries(_NAMES, st.tuples(ports, ports), min_size=1, max_size=3))
    found = {stem: (set(ins), set(outs)) for stem, (ins, outs) in drawn.items()}
    for x, separator, y, z in draw(st.lists(st.tuples(_NAMES, _SEPARATORS, _NAMES, _NAMES), max_size=2)):
        for stem, port in ((x, f'{y}{separator}{z}'), (f'{x}{separator}{y}', z)):
            for side in found.setdefault(stem, (set(), set())):
                side.add(port)
    return {stem: (frozenset(ins), frozenset(outs)) for stem, (ins, outs) in found.items()}


@st.composite
def _trees(draw: st.DrawFn) -> tuple[str, Yaml, dict[str, _Signature]]:
    """A root name, a workflow, and the signatures of the tools it runs.

    The root runs every tool at least once, so names built to meet do meet,
    and may repeat a step `x` beside a step named `x·n`, `n` being the
    position of one copy of `x`. Steps repeat, a subworkflow may be called
    more than once, and a call may be named like a tool; subworkflows nest up
    to two levels below the root.
    A subworkflow name keeps one body throughout, since the registry holds
    one body per name.
    """
    signatures = draw(_signatures())
    stems = st.lists(st.sampled_from(sorted(signatures)), min_size=1, max_size=3)
    callees = draw(st.lists(_NAMES, max_size=2, unique=True))
    bodies: dict[str, Yaml] = {}

    def calls(names: list[str]) -> list[str]:
        return draw(st.lists(st.sampled_from(names), max_size=2)) if names else []

    for depth, name in reversed(list(enumerate(callees))):
        bodies[name] = _workflow(draw(stems), signatures, calls(callees[depth + 1:]), bodies)
    order = draw(st.permutations(sorted(signatures)))
    repeated = draw(st.lists(st.tuples(st.sampled_from(order), _SEPARATORS), max_size=2, unique_by=lambda x: x[0]))
    numbered = sorted({f'{x}{separator}{order.index(x) + 1}' for x, separator in repeated} - set(signatures))
    signatures |= {name: (frozenset(), frozenset()) for name in numbered}
    root = order + numbered + [x for x, _separator in repeated] + draw(stems)
    return draw(_NAMES), _workflow(root, signatures, calls(callees), bodies), signatures


def _exposed(yml: Yaml, signatures: dict[str, _Signature]) -> tuple[int, int]:
    """How many inputs and outputs `yml` exposes: one per literal and one per
    tool output, through every subworkflow call. Counted from the source, so
    a made name that two of them share shows as a missing key."""
    inputs = outputs = 0
    for step in yml['steps']:
        if 'subtree' in step:
            nested = _exposed(step['subtree'], signatures)
            inputs, outputs = inputs + nested[0], outputs + nested[1]
        else:
            inputs, outputs = inputs + len(step['in']), outputs + len(signatures[step['id']][1])
    return inputs, outputs


def _assert_every_name_is_distinct(yml: Yaml, signatures: dict[str, _Signature],
                                   artifact: CompilationArtifact, files: list[str]) -> None:
    """Each step of `yml` keeps its own id, each exposed port its own key, and
    each process its own file, in `artifact` and every workflow it calls."""
    cwl = artifact.cwl
    ids = [str(step['id']) for step in cwl['steps']]
    assert len(set(ids)) == len(yml['steps']) == len(ids), f'steps share an id in {artifact.name}: {ids}'
    inputs, outputs = _exposed(yml, signatures)
    assert len(cwl['inputs']) == inputs, f'inputs share a name in {artifact.name}: {sorted(cwl["inputs"])}'
    assert len(cwl['outputs']) == outputs, f'outputs share a name in {artifact.name}: {sorted(cwl["outputs"])}'
    assert len(artifact.children) == len(yml['steps'])
    for step, child in zip(yml['steps'], artifact.children):
        files.append(child.name)
        if 'subtree' in step:
            _assert_every_name_is_distinct(step['subtree'], signatures, child, files)


@pytest.mark.skip_pypi_ci
@pytest.mark.slow
@given(_trees())
@example(('r', {'steps': [{'id': 'a', 'in': {'b_c': {'wic_inline_input': 'v'}}},
                          {'id': 'a_b', 'in': {'c': {'wic_inline_input': 'w'}}}]},
          {'a': (frozenset({'b_c'}), frozenset({'b_c'})), 'a_b': (frozenset({'c'}), frozenset({'c'}))}))
@ORACLE
def test_distinct_authored_names_never_share_a_made_name(tree: tuple[str, Yaml, dict[str, _Signature]]) -> None:
    """No two distinct steps or ports get the same made name, whatever
    `[A-Za-z0-9_.-]` names they were written with.

    A made name is a step id, `<step>~<port>` for a port exposed through a
    step, or `<root>~<address>` for a written file; a step id repeated in one
    workflow is `<id>@<index>` on every copy. The drawn names are built to
    meet under a join on `_`, `.` or `-`, steps repeat, a call may share a
    tool's name, and subworkflows nest and repeat, so each of those spellings
    is exercised where it could collide. A collision is a dict key silently
    overwritten, which the compiled document cannot show by itself; the
    expected counts are therefore taken from the source.
    """
    root, yml, signatures = tree
    compiled = compile_hermetic(copy.deepcopy(yml), root, tools=_tools(signatures))
    files: list[str] = []
    _assert_every_name_is_distinct(yml, signatures, compiled.artifact, files)
    duplicates = sorted({name for name in files if files.count(name) > 1})
    assert not duplicates, f'two processes are written to the same file: {duplicates}'
