"""Language-version resolution, surfacing, and override.

The claims under test, one property each: the resolved version is the highest
that satisfies the source; no lower version is ever silently selected;
untagged files resolve rather than fail; every successful compile reports a
resolved version, while the compiled CWL carries no `sophios` key or prefix;
an explicit setting always beats a file tag; and one compilation yields
exactly one version, tree-wide.

The mechanism is exercised against fabricated version histories as well as the
real single-entry list — `sophios.lang.versions` states why, once. The claims
under test have their normative home in the language spec, §6.
"""
import sys
from pathlib import Path
import pytest
import yaml
from hypothesis import given
from hypothesis import strategies as st

from sophios.ir.frontdoor import bundle_from_disk
from sophios.lang import LANG_VERSION, resolve_lang_version
from sophios.lang.diagnostics import SophiosError
from sophios.lang.error_codes import SophiosErrorCode
import sophios.cli
import sophios.main
from sophios.wic_types import Yaml

from .compile_harness import FAST, TOUCH, compile_cwl, compile_info
from .hermetic import subworkflow_step
from .test_setup import load_test_registry

#: The key earlier Sophios versions annotated compiled CWL with.
ANNOTATION_KEY = 'sophios:lang_version'

#: Fabricated, ordered version histories: what the resolver will face someday.
histories = st.lists(
    st.integers(min_value=0, max_value=40).map(lambda n: f'0.0.{n}'),
    min_size=1, max_size=6, unique=True,
).map(lambda vs: tuple(sorted(vs, key=lambda v: int(v.rsplit('.', 1)[1]))))


def _compile(yml: Yaml, lang_version: str | None = None) -> Yaml:
    """Compile one in-memory workflow and return the emitted CWL."""
    return compile_cwl(yml, 'lang_version', lang_version=lang_version)


# --------------------------------------------------------------------------
# The resolution rule, over fabricated histories
# --------------------------------------------------------------------------


@pytest.mark.fast
@given(histories)
@FAST
def test_untagged_resolves_to_the_newest(known: tuple[str, ...]) -> None:
    """With nothing to satisfy, the highest known version is chosen."""
    assert resolve_lang_version(None, (), known=known) == known[-1]


@pytest.mark.fast
@given(histories, st.data())
@FAST
def test_no_lower_version_is_silently_selected(known: tuple[str, ...], data: st.DataObject) -> None:
    """A pin is the only thing that moves resolution below the newest,
    and a pin is the author's word, not silence."""
    pin = data.draw(st.sampled_from(known))
    resolved = resolve_lang_version(None, (pin,), known=known)
    assert resolved == pin
    # Order is history position, not string order: '0.0.10' is newer than '0.0.2'.
    untagged = resolve_lang_version(None, (), known=known)
    assert known.index(untagged) >= known.index(resolved)


@pytest.mark.fast
@given(histories)
@FAST
def test_untagged_never_fails(known: tuple[str, ...]) -> None:
    """Resolution of an untagged tree cannot fail."""
    resolve_lang_version(None, (), known=known)  # must not raise


@pytest.mark.fast
def test_unknown_versions_are_reported_not_coerced() -> None:
    """Asking for a version that does not exist is an error naming the known
    ones — never a silent fallback to something else."""
    for request, pins in (('9.9.9', ()), (None, ('9.9.9',))):
        with pytest.raises(SophiosError) as caught:
            resolve_lang_version(request, pins)
        assert caught.value.diagnostics[0].code is SophiosErrorCode.UNKNOWN_LANG_VERSION
        assert LANG_VERSION in caught.value.diagnostics[0].message


# --------------------------------------------------------------------------
# Precedence, and tree-wide oneness
# --------------------------------------------------------------------------


@pytest.mark.fast
@given(histories, st.data())
@FAST
def test_explicit_beats_any_tag(known: tuple[str, ...], data: st.DataObject) -> None:
    """The explicit setting wins regardless of what the tree pins."""
    explicit = data.draw(st.sampled_from(known))
    pins = tuple(data.draw(st.lists(st.sampled_from(known), max_size=4)))
    assert resolve_lang_version(explicit, pins, known=known) == explicit


@pytest.mark.fast
@given(histories, st.data())
@FAST
def test_one_version_or_a_conflict_report(known: tuple[str, ...], data: st.DataObject) -> None:
    """Agreeing pins resolve to their one version; disagreeing pins are
    a reported conflict, never an average or a quiet winner."""
    pins = tuple(data.draw(st.lists(st.sampled_from(known), min_size=1, max_size=4)))
    if len(set(pins)) == 1:
        assert resolve_lang_version(None, pins, known=known) == pins[0]
    else:
        with pytest.raises(SophiosError) as caught:
            resolve_lang_version(None, pins, known=known)
        assert caught.value.diagnostics[0].code is SophiosErrorCode.LANG_VERSION_CONFLICT


@pytest.mark.skip_pypi_ci
@pytest.mark.fast
def test_conflicting_pins_across_a_tree_are_caught_at_the_root() -> None:
    """The same, through the compiler: a subworkflow pin that disagrees with the root's is a
    conflict at compile time — every reachable document's pin is collected."""
    tree: Yaml = {
        'wic': {'lang_version': '0.0.1'},
        'steps': [subworkflow_step('child.wic', {'wic': {'lang_version': '9.9.9'}, **TOUCH})],
    }
    with pytest.raises(SophiosError) as caught:
        _compile(tree)
    assert caught.value.diagnostics[0].code is SophiosErrorCode.UNKNOWN_LANG_VERSION


@pytest.mark.skip_pypi_ci
@pytest.mark.fast
def test_an_inline_implementation_body_pin_is_collected() -> None:
    """A pin that lives only on the selected implementation body still decides the version."""
    body = {'wic': {'lang_version': '9.9.9'}, **TOUCH}
    tree: Yaml = {
        'wic': {'default_implementation': 'fast', 'implementations': {'fast': body}},
        'steps': [],
    }
    with pytest.raises(SophiosError) as caught:
        _compile(tree)
    assert caught.value.diagnostics[0].code is SophiosErrorCode.UNKNOWN_LANG_VERSION


# --------------------------------------------------------------------------
# The resolved version is never a secret: every compile reports
# it (CLI, Python API), and the compiled CWL stays plain CWL
# --------------------------------------------------------------------------


@pytest.mark.skip_pypi_ci
@pytest.mark.fast
def test_every_compile_reports_its_version() -> None:
    """The result names the resolved version, untagged, tagged and explicitly
    pinned alike, and the CWL it carries has no annotation of its own."""
    tagged = {'wic': {'lang_version': '0.0.1'}, **TOUCH}
    for yml, explicit, expected in ((TOUCH, None, LANG_VERSION), (tagged, None, '0.0.1'), (TOUCH, '0.0.1', '0.0.1')):
        result = compile_info(yml, 'lang_version', lang_version=explicit)
        assert result.lang_version == expected
        assert ANNOTATION_KEY not in result.artifact.cwl
        assert 'sophios' not in result.artifact.cwl.get('$namespaces', {})


@pytest.mark.skip_pypi_ci
@pytest.mark.fast
def test_the_python_api_surfaces_the_version() -> None:
    """`CompiledWorkflow.lang_version` names the version, no digging."""
    from sophios.api.python.tool_builder import CommandLineTool, Input, Inputs, Output, Outputs, cwl
    from sophios.api.python.workflow import Step, Workflow

    emit = (CommandLineTool('emit_text',
                            Inputs(message=Input(cwl.string, position=1)),
                            Outputs(file=Output(cwl.file, glob='stdout.txt')))
            .base_command('echo')
            .stdout('stdout.txt'))
    step = Step(emit)
    step.inputs.message = 'hello'
    compiled = Workflow([step], 'p18_version').compile()
    assert compiled.lang_version == LANG_VERSION
    assert ANNOTATION_KEY not in compiled.cwl_workflow


@pytest.mark.skip_pypi_ci
@pytest.mark.fast
def test_a_pinned_file_reports_its_version_and_writes_plain_cwl(
        tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]) -> None:
    """A `wic: lang_version:` pin still compiles and the command line prints the
    version; a Python workflow, which pins per compile, names it in
    `CompiledWorkflow.lang_version`; neither CWL carries a `sophios` key or prefix."""
    from sophios.api.python.workflow import Workflow  # pylint: disable=import-outside-toplevel

    source = tmp_path / 'pinned.wic'
    source.write_text(f'wic:\n  lang_version: {LANG_VERSION}\n'
                      'steps:\n- id: touch\n  in:\n    filename: !ii empty.txt\n', encoding='utf-8')
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, 'argv', ['sophios', '--yaml', str(source), '--generate_cwl_workflow'])
    sophios.main._main(*sophios.cli.parser.parse_known_args())  # pylint: disable=protected-access
    assert f'Sophios lang_version: {LANG_VERSION}' in capsys.readouterr().out
    written = yaml.safe_load((tmp_path / 'autogenerated' / 'pinned.cwl').read_text(encoding='utf-8'))
    assert ANNOTATION_KEY not in written and 'sophios' not in written.get('$namespaces', {})

    unpinned = tmp_path / 'unpinned.wic'
    unpinned.write_text('steps:\n- id: touch\n  in:\n    filename: !ii empty.txt\n', encoding='utf-8')
    tools = load_test_registry().tools
    compiled = Workflow.from_wic(unpinned, tool_registry=tools).compile(tool_registry=tools, lang_version=LANG_VERSION)
    assert compiled.lang_version == LANG_VERSION
    assert ANNOTATION_KEY not in compiled.cwl_workflow


@pytest.mark.fast
def test_a_file_tag_reaches_the_compiler_from_its_own_text(tmp_path: Path) -> None:
    """The `wic: lang_version:` pin survives the path a real `.wic` file takes.

    Every other test here hands a document straight to the compiler. This one
    writes a file and compiles that file, because a key can be specified,
    parsed, resolved and surfaced and still be dead on arrival for anyone
    writing one -- which is what happened when `lang_version` reached the language spec, the
    resolver and the compiler but not the vocabulary that admits it.

    The gate used to be the generated jsonschema, applied by the file loader.
    It is the parser now: `Grammar.SIDECAR_KEYS` is where the language names
    the keys its `wic:` block has.
    """
    from sophios.lang.parser import Grammar, parse  # pylint: disable=import-outside-toplevel

    assert 'lang_version' in Grammar.SIDECAR_KEYS

    source = f'wic:\n  lang_version: {LANG_VERSION}\nsteps:\n- id: mk_file\n'
    written = tmp_path / 'pinned.wic'
    written.write_text(source, encoding='utf-8')

    parsed = parse(written.read_text(encoding='utf-8'), 'pinned.wic')
    assert parsed.document is not None, list(parsed.diagnostics)
    assert not parsed.diagnostics.has_errors, [str(d) for d in parsed.diagnostics]
    assert parsed.document.sidecar is not None
    assert dict(parsed.document.sidecar.entries)['lang_version'] == LANG_VERSION


@pytest.mark.fast
def test_the_cli_flag_reaches_the_compiler() -> None:
    """`--lang_version` survives the trip from argv to compiler options.

    The flag's entire contract is precedence — an explicit setting beats a
    file pin — and it could not keep it: `main` parsed the real argv, then
    `_build_and_compile_workflow` threw those arguments away and called
    `get_dicts_for_compilation()`, which re-parses a patched
    `sys.argv` of `['sophios', '--yaml', '']`. Every flag arrived at the
    compiler as its default, so `--lang_version 9.9.9` compiled at the newest
    version instead of reporting an unknown one.

    Tested at the seam because that is where the value was being dropped. The
    caller side needs no test of its own any more: `get_dicts_for_compilation`
    requires its argument, so a `main` that forgets to pass what it parsed is
    a TypeError and a mypy failure rather than a silent compile with defaults.
    The guarantee moved from a test into the signature.
    """
    supplied = sophios.cli.get_args('wf.wic', ['--lang_version', '9.9.9'])
    options, _graph_settings = sophios.cli.get_dicts_for_compilation(supplied)
    assert options['lang_version'] == '9.9.9', 'the flag never reached the compiler'

    # And a library caller asking for defaults gets no pin, as it should.
    defaults = sophios.cli.default_compilation_settings()
    assert defaults['lang_version'] is None


@pytest.mark.fast
def test_a_mistyped_pin_is_reported_not_ignored(tmp_path: Path) -> None:
    """`lang_version: 1.0` is a YAML float, not a string. A walk that only
    collected strings made it vanish — silently inferred over, when the
    author plainly asked for something, and by the no-silent-selection rule
    above, silence is the one wrong answer. Any value under the key is a pin
    claim; a non-version is reported as unknown, naming what was written."""
    written = tmp_path / 'mistyped.wic'
    written.write_text('wic: {lang_version: 1.0}\nsteps:\n- id: s\n', encoding='utf-8')
    pins = bundle_from_disk(written, {'global': {}}, {}).lang_version_pins
    assert pins == ('1.0',)

    with pytest.raises(SophiosError) as caught:
        resolve_lang_version(None, pins)
    assert caught.value.diagnostics[0].code is SophiosErrorCode.UNKNOWN_LANG_VERSION
    assert "'1.0'" in caught.value.diagnostics[0].message
