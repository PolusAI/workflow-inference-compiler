"""A tool that refers to files beside it is written packed, so it means the same wherever it is written."""
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any, Final

import pytest
import yaml
from schema_salad.fetcher import DefaultFetcher

import sophios
from sophios import main, run_local, utils_cwl
from sophios.api.python import _workflow_runtime as python_runtime
from sophios.api.python.workflow import Step, Workflow
from sophios.cli import default_compilation_settings
from sophios.compiler import compile_source
from sophios.ir.artifacts import CompilationArtifact
from sophios.ir.frontdoor import bundle_from_disk
from sophios.lang.diagnostics import SophiosError
from sophios.lang.error_codes import SophiosErrorCode
from sophios.wic_types import StepId, Tool, Tools

#: Its script by `$include`, its data by a relative default: both beside it.
RELCAT: Final = '''cwlVersion: v1.2
class: CommandLineTool
requirements:
  InitialWorkDirRequirement:
    listing:
    - entryname: s.sh
      entry:
        $include: s.sh
baseCommand: [sh, s.sh]
inputs:
  extra:
    type: File
    default: {class: File, location: data.txt}
    inputBinding: {position: 1}
outputs:
  out: {type: stdout}
stdout: out.txt
'''
SCRIPT: Final = 'echo from-include; cat "$1"\n'
DATA: Final = 'data beside the tool\n'
#: What `relcat` prints.
PRINTED: Final = 'from-include\n' + DATA
PLAIN: Final = '''cwlVersion: v1.2
class: CommandLineTool
baseCommand: echo
inputs: {message: {type: string, default: hi, inputBinding: {position: 1}}}
outputs: {out: {type: stdout}}
stdout: out.txt
'''
#: A relative default, and an output whose format inference reads as written: `edam:format_2330`.
MAKE: Final = '''cwlVersion: v1.2
class: CommandLineTool
$namespaces: {edam: https://edamontology.org/}
baseCommand: cat
stdout: made.txt
inputs:
  seed: {type: File, default: {class: File, location: seed.txt}, inputBinding: {position: 1}}
outputs:
  out: {type: File, format: edam:format_2330, outputBinding: {glob: made.txt}}
'''
USE: Final = '''cwlVersion: v1.2
class: CommandLineTool
$namespaces: {edam: https://edamontology.org/}
baseCommand: cat
stdout: used.txt
inputs:
  file: {type: File, format: edam:format_2330, inputBinding: {position: 1}}
outputs:
  out: {type: File, outputBinding: {glob: used.txt}}
'''
MAKE_USE: Final = 'steps:\n- id: make\n  run: tools/make.cwl\n- id: use\n  run: tools/use.cwl\n'
#: Like `extract_pdbbind_refined.cwl`: its image is built from a Dockerfile beside it.
DOCKERFILE_TOOL: Final = '''cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerImageId: built_here
    dockerFile:
      $include: Dockerfile_built
baseCommand: echo
inputs: {}
outputs: {}
'''
DOCKERFILE: Final = 'FROM alpine:3.20\nRUN echo built\n'
#: Its `$include` names nothing.
BROKEN: Final = '''cwlVersion: v1.2
class: CommandLineTool
requirements:
  InitialWorkDirRequirement:
    listing: [{entryname: s.sh, entry: {$include: missing.sh}}]
baseCommand: [sh, s.sh]
inputs: {}
outputs: {}
'''


def _relcat(directory: Path, script: str = SCRIPT) -> Path:
    """`relcat.cwl`, its script and its data, in `directory`."""
    directory.mkdir(parents=True, exist_ok=True)
    (directory / 's.sh').write_text(script, encoding='utf-8')
    (directory / 'data.txt').write_text(DATA, encoding='utf-8')
    tool = directory / 'relcat.cwl'
    tool.write_text(RELCAT, encoding='utf-8')
    return tool


def _registry(*tools: Path) -> Tools:
    """Each file as a search-path tool, read as the loader reads it."""
    return {StepId(tool.stem, 'global'): Tool(str(tool), yaml.safe_load(tool.read_text(encoding='utf-8')))
            for tool in tools}


def _compile(wic: Path, tools: Tools | None = None) -> CompilationArtifact:
    return compile_source(bundle_from_disk(wic, {}, tools or {}), default_compilation_settings(),
                          testing=True).artifact


def _tool(artifact: CompilationArtifact, step: str) -> dict[str, Any]:
    """The description the artifact carries for `step`'s tool."""
    child, = (child for child in artifact.children if child.namespace == (step,))
    return child.cwl


def _requirement(tool: dict[str, Any], section: str, name: str) -> Any:
    """A requirement or hint by its class, whether the section is written as a mapping or as a list."""
    entries = tool[section]
    if isinstance(entries, list):
        entries = {item['class']: {key: value for key, value in item.items() if key != 'class'} for item in entries}
    return entries[name]


def _entry(tool: dict[str, Any]) -> Any:
    """`relcat`'s script entry."""
    return _requirement(tool, 'requirements', 'InitialWorkDirRequirement')['listing'][0]['entry']


@pytest.mark.fast
def test_a_tool_with_an_include_and_a_relative_default_is_emitted_packed(tmp_path: Path) -> None:
    """The include is its text and the default an absolute URI; the registry still holds the tool as read."""
    tool = _relcat(tmp_path / 'tools')
    (tmp_path / 'w.wic').write_text('steps:\n- id: relcat\n', encoding='utf-8')
    registry = _registry(tool)

    written = _tool(_compile(tmp_path / 'w.wic', registry), 'relcat')

    assert _entry(written) == SCRIPT
    assert written['inputs']['extra']['default'] == {
        'class': 'File', 'location': (tmp_path / 'tools/data.txt').as_uri()}
    assert written['cwlVersion'] == 'v1.2'
    assert _entry(registry[StepId('relcat', 'global')].cwl) == {'$include': 's.sh'}


@pytest.mark.fast
def test_a_tool_without_relative_references_is_emitted_as_read(tmp_path: Path) -> None:
    """Nothing to resolve, nothing packed: the description is the tool as the loader read it."""
    (tmp_path / 'echo.cwl').write_text(PLAIN, encoding='utf-8')
    (tmp_path / 'w.wic').write_text('steps:\n- id: echo\n', encoding='utf-8')
    registry = _registry(tmp_path / 'echo.cwl')

    assert _tool(_compile(tmp_path / 'w.wic', registry), 'echo') == registry[StepId('echo', 'global')].cwl


@pytest.mark.fast
def test_a_tool_three_steps_run_is_packed_once(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """One compile packs a tool once, however many steps run it, and every step's description is packed."""
    tool = _relcat(tmp_path / 'tools')
    (tmp_path / 'w.wic').write_text('steps:\n- id: relcat\n- id: relcat\n- id: relcat\n', encoding='utf-8')
    packs: list[str] = []
    pack = utils_cwl.packed

    def counted(text: str, base_uri: str) -> dict[str, Any]:
        packs.append(base_uri)
        return pack(text, base_uri)

    monkeypatch.setattr(utils_cwl, 'packed', counted)

    artifact = _compile(tmp_path / 'w.wic', _registry(tool))

    assert packs == [tool.as_uri()]
    assert [_entry(child.cwl) for child in artifact.children] == [SCRIPT] * 3


@pytest.mark.fast
def test_an_inline_body_with_a_relative_reference_is_packed_against_its_wic(tmp_path: Path,
                                                                            monkeypatch: pytest.MonkeyPatch) -> None:
    """An inline `run:` body is written in its `.wic`, so `s.sh` is the one beside the `.wic`, not the working
    directory's."""
    project = tmp_path / 'project'
    project.mkdir()
    (project / 's.sh').write_text(SCRIPT, encoding='utf-8')
    (project / 'w.wic').write_text('''steps:
- id: inline
  run:
    class: CommandLineTool
    requirements:
      InitialWorkDirRequirement:
        listing: [{entryname: s.sh, entry: {$include: s.sh}}]
    baseCommand: [sh, s.sh]
    inputs: {}
    outputs: {}
''', encoding='utf-8')
    (tmp_path / 'elsewhere').mkdir()
    (tmp_path / 'elsewhere' / 's.sh').write_text('echo the wrong one\n', encoding='utf-8')
    monkeypatch.chdir(tmp_path / 'elsewhere')

    assert _entry(_tool(_compile(project / 'w.wic'), 'inline')) == SCRIPT


@pytest.mark.fast
def test_one_inline_body_in_two_directories_reads_each_ones_files(tmp_path: Path) -> None:
    """The same body, written in two `.wic` files in two directories, includes each directory's `s.sh`."""
    body = '''- id: say
  run:
    class: CommandLineTool
    requirements:
      InitialWorkDirRequirement:
        listing: [{entryname: s.sh, entry: {$include: s.sh}}]
    baseCommand: [sh, s.sh]
    inputs: {}
    outputs: {}
'''
    (tmp_path / 'a').mkdir()
    (tmp_path / 'b').mkdir()
    (tmp_path / 'a' / 's.sh').write_text('echo a\n', encoding='utf-8')
    (tmp_path / 'b' / 's.sh').write_text('echo b\n', encoding='utf-8')
    (tmp_path / 'b' / 'inner.wic').write_text('steps:\n' + body, encoding='utf-8')
    (tmp_path / 'a' / 'w.wic').write_text('steps:\n' + body + '- id: inner\n  run: ../b/inner.wic\n',
                                          encoding='utf-8')

    artifact = _compile(tmp_path / 'a' / 'w.wic')

    inner, = (child for child in artifact.children if child.namespace == ('inner',))
    assert (_entry(_tool(artifact, 'say')), _entry(_tool(inner, 'say'))) == ('echo a\n', 'echo b\n')


def _dockerfile_project(tmp_path: Path) -> Path:
    """A project whose search path holds `built.cwl` and the Dockerfile it builds its image from."""
    tools = tmp_path / 'tools'
    tools.mkdir()
    (tools / 'built.cwl').write_text(DOCKERFILE_TOOL, encoding='utf-8')
    (tools / 'Dockerfile_built').write_text(DOCKERFILE, encoding='utf-8')
    (tmp_path / 'w.wic').write_text('steps:\n- id: built\n', encoding='utf-8')
    (tmp_path / 'config.json').write_text(json.dumps({'search_paths_cwl': {'global': [str(tools)]},
                                                      'search_paths_wic': {'global': [str(tmp_path)]}}),
                                          encoding='utf-8')
    return tmp_path


@pytest.mark.fast
def test_a_dockerfile_include_is_packed_text(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """The written tool carries the Dockerfile itself, in the tree and embedded, and keeps its own `v1.0`."""
    monkeypatch.chdir(_dockerfile_project(tmp_path))
    for flags in ((), ('--cwl_inline_runtag',)):
        monkeypatch.setattr('sys.argv', ['sophios', '--yaml', 'w.wic', '--config_file', 'config.json',
                                         '--generate_cwl_workflow', '--basepath', 'out', *flags])
        main.main()
        root = yaml.safe_load(Path('out/w.cwl').read_text(encoding='utf-8'))
        tool = root['steps'][0]['run'] if flags else yaml.safe_load(Path('out/w~built.cwl').read_text(encoding='utf-8'))
        docker = _requirement(tool, 'hints', 'DockerRequirement')
        assert docker == {'dockerImageId': 'built_here', 'dockerFile': DOCKERFILE}, flags
        assert flags or tool['cwlVersion'] == 'v1.0'


def _analysis_project(tmp_path: Path) -> Path:
    """`relcat` on the search path, run by `w.wic` and declared as a real-time analysis by `live.wic`."""
    tools = tmp_path / 'tools'
    _relcat(tools)
    (tools / 'slow_log.cwl').write_text('''cwlVersion: v1.2
class: CommandLineTool
baseCommand: [touch, out.log]
inputs: {}
outputs: {log: {type: File, outputBinding: {glob: out.log}}}
''', encoding='utf-8')
    adapters = Path(__file__).resolve().parents[2] / 'cwl_adapters' / 'file_watchers'
    (tmp_path / 'w.wic').write_text('steps:\n- id: relcat\n', encoding='utf-8')
    (tmp_path / 'live.wic').write_text('''steps:
  - id: slow_log
  - id: cwl_subinterpreter
    in:
      file_pattern: !ii out.log
      cwl_tool: !ii relcat
      max_times: !ii '1'
      config: !ii {in: {}}
''', encoding='utf-8')
    (tmp_path / 'config.json').write_text(json.dumps({'search_paths_cwl': {'global': [str(tools), str(adapters)]},
                                                      'search_paths_wic': {'global': [str(tmp_path)]}}),
                                          encoding='utf-8')
    return tools / 'relcat.cwl'


@pytest.mark.fast
def test_compile_run_and_an_analysis_all_carry_the_packed_text(tmp_path: Path,
                                                               monkeypatch: pytest.MonkeyPatch) -> None:
    """Python `compile()`, the document `run()` writes, and a real-time analysis's tool are packed alike."""
    tool = _analysis_project(tmp_path)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(run_local, 'generate_run_script', lambda _cmdline: None)
    monkeypatch.setattr(python_runtime.preflight, 'prepare', lambda *_a, **_k: None)

    compiled = Workflow([Step(clt_path=tool)], 'py').compile().cwl_workflow
    Workflow([Step(clt_path=tool)], 'py').run(basepath='ran', run_args_dict={'generate_run_script': 'yes'})
    ran = yaml.safe_load(Path('ran/py.cwl').read_text(encoding='utf-8'))
    monkeypatch.setattr('sys.argv', ['sophios', '--yaml', 'live.wic', '--config_file', 'config.json',
                                     '--generate_cwl_workflow'])
    main.main()
    analysis, = Path('autogenerated/realtime/live').glob('*~relcat.cwl')

    for where, written in (('compile()', compiled['steps'][0]['run']), ('run()', ran['steps'][0]['run']),
                           ('analysis', yaml.safe_load(analysis.read_text(encoding='utf-8')))):
        assert _entry(written) == SCRIPT, where
        assert written['inputs']['extra']['default']['location'] == (tool.parent / 'data.txt').as_uri(), where


@pytest.mark.fast
def test_inference_wires_a_packed_tool_as_it_wires_one_as_written(tmp_path: Path) -> None:
    """Inference reads `make` as written, so its `edam:format_2330` output feeds `use`'s input: no new root
    input. The description written for `make` is packed."""
    tools = tmp_path / 'tools'
    tools.mkdir()
    (tools / 'seed.txt').write_text('seed\n', encoding='utf-8')
    (tools / 'make.cwl').write_text(MAKE, encoding='utf-8')
    (tools / 'use.cwl').write_text(USE, encoding='utf-8')
    (tmp_path / 'w.wic').write_text(MAKE_USE, encoding='utf-8')

    artifact = _compile(tmp_path / 'w.wic')

    use, = (step for step in artifact.cwl['steps'] if step['id'] == 'use')
    assert use['in'] == {'file': 'make/out'}
    assert not artifact.cwl['inputs']
    assert _tool(artifact, 'make')['inputs']['seed']['default']['location'] == (tools / 'seed.txt').as_uri()


@pytest.mark.fast
def test_a_broken_tool_nobody_runs_stops_no_compile(tmp_path: Path) -> None:
    """Only a tool a step runs is packed: one on the search path that cannot be packed costs nothing."""
    (tmp_path / 'broken.cwl').write_text(BROKEN, encoding='utf-8')
    (tmp_path / 'echo.cwl').write_text(PLAIN, encoding='utf-8')
    (tmp_path / 'w.wic').write_text('steps:\n- id: echo\n', encoding='utf-8')

    artifact = _compile(tmp_path / 'w.wic', _registry(tmp_path / 'broken.cwl', tmp_path / 'echo.cwl'))

    assert [child.namespace for child in artifact.children] == [('echo',)]


@pytest.mark.fast
def test_a_step_running_a_tool_with_a_missing_include_reports_cwltools_error_at_the_step(tmp_path: Path) -> None:
    """cwltool's own message, at the line of the step that runs the tool."""
    (tmp_path / 'broken.cwl').write_text(BROKEN, encoding='utf-8')
    (tmp_path / 'echo.cwl').write_text(PLAIN, encoding='utf-8')
    (tmp_path / 'w.wic').write_text('steps:\n- id: echo\n- id: broken\n', encoding='utf-8')

    with pytest.raises(SophiosError) as raised:
        _compile(tmp_path / 'w.wic', _registry(tmp_path / 'broken.cwl', tmp_path / 'echo.cwl'))

    diagnostic, = raised.value.diagnostics
    assert diagnostic.code is SophiosErrorCode.SUBWORKFLOW_INVALID
    assert diagnostic.locator is not None and (diagnostic.locator.step, diagnostic.locator.index) == ('broken', 2)
    assert diagnostic.span is not None and diagnostic.span.start_line == 3
    assert str(tmp_path / 'broken.cwl') in diagnostic.message
    # The OSError as schema-salad passes it on. On Windows the error's filename, spelled with backslashes,
    # differs from the URI's path, so schema-salad prefixes it with `Error reading <uri>: `.
    reason = diagnostic.message.splitlines()[-1].strip()
    assert "[Errno 2] No such file or directory: '" in reason and reason.endswith("missing.sh'")


def _run(workflow: Path, *flags: str) -> None:
    """`sophios --run_local` on `workflow` from the current directory, in its own process as on the command
    line: Toil and cwltool share logging state that one process running both leaves broken."""
    package = str(Path(sophios.__file__).resolve().parents[1])
    env = {**os.environ, 'PYTHONPATH': os.pathsep.join(filter(None, (package, os.environ.get('PYTHONPATH'))))}
    subprocess.run([sys.executable, '-m', 'sophios.main', '--yaml', str(workflow), '--run_local',
                    '--ignore_docker_install', '--ignore_docker_processes', *flags], check=True, env=env)


def _produced() -> list[str]:
    """What the run's output files hold, in order: cwltool and Toil each write them to an outdir of their own."""
    return sorted(path.read_text(encoding='utf-8') for path in Path('autogenerated').glob('outdir_*/**/*')
                  if path.is_file())


RUNNERS: Final = ['cwltool', 'toil-cwl-runner']
FORMS: Final = {'tree': (), 'embedded': ('--cwl_inline_runtag',)}


@pytest.mark.needs_cwltool
@pytest.mark.parametrize('runner', RUNNERS)
@pytest.mark.parametrize('form', FORMS)
def test_a_packed_tool_runs_from_the_tree_and_embedded_twice(tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
                                                             runner: str, form: str) -> None:
    """Two steps run `relcat`, from the tree's `w~<step>.cwl` and from one embedded document: each prints its
    included script's line and its default's data."""
    _relcat(tmp_path / 'tools')
    (tmp_path / 'w.wic').write_text('steps:\n- id: one\n  run: tools/relcat.cwl\n'
                                    '- id: two\n  run: tools/relcat.cwl\n', encoding='utf-8')
    monkeypatch.chdir(tmp_path)

    _run(tmp_path / 'w.wic', '--cwl_runner', runner, *FORMS[form])

    assert _produced() == [PRINTED, PRINTED]


SCHEMADEF: Final = Path(__file__).resolve().parent / 'data' / 'cwl-v1.2' / 'tests'


@pytest.mark.needs_cwltool
def test_a_schemadef_tool_with_an_import_is_written_runnable(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """The conformance tool whose record type is `$import`ed from a file beside it: the description the tree
    holds for its step runs, given the conformance job."""
    tool = SCHEMADEF / 'schemadef_types_with_import-tool.cwl'
    job = json.loads((SCHEMADEF / 'schemadef_types_with_import-job.json').read_text(encoding='utf-8'))
    value = job['readgroups_bam']
    value['bam'] = {'class': 'File', 'location': (SCHEMADEF / value['bam']['location']).as_uri()}
    (tmp_path / 'job.yml').write_text(yaml.safe_dump({'message': value}), encoding='utf-8')
    (tmp_path / 'w.wic').write_text(yaml.safe_dump({'steps': [{'id': 'sd', 'run': os.path.relpath(tool, tmp_path)}]}),
                                    encoding='utf-8')
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr('sys.argv', ['sophios', '--yaml', 'w.wic', '--generate_cwl_workflow'])
    main.main()

    ran = subprocess.run([sys.executable, '-c', 'import sys, cwltool.main; sys.exit(cwltool.main.main())',
                          '--outdir', 'out', 'autogenerated/w~sd.cwl', 'job.yml'], capture_output=True, text=True,
                         check=False)

    assert ran.returncode == 0, ran.stderr
    assert json.loads(ran.stdout)['out'] == '1970-01-01T00:00:00Z\n'


@pytest.mark.needs_cwltool
def test_a_packed_tool_keeps_the_partial_failure_rewrite_written_and_embedded(tmp_path: Path,
                                                                              monkeypatch: pytest.MonkeyPatch) -> None:
    """`--partial_failure_enable` rewrites the packed description: its output optional, exit 1 a success, in the
    tree and embedded; the tree runs."""
    _relcat(tmp_path / 'tools')
    (tmp_path / 'w.wic').write_text('steps:\n- id: relcat\n  run: tools/relcat.cwl\n', encoding='utf-8')
    monkeypatch.chdir(tmp_path)

    _run(tmp_path / 'w.wic', '--partial_failure_enable')
    tree = yaml.safe_load(Path('autogenerated/w~relcat.cwl').read_text(encoding='utf-8'))
    assert _produced() == [PRINTED]
    _run(tmp_path / 'w.wic', '--partial_failure_enable', '--cwl_inline_runtag')
    embedded = yaml.safe_load(Path('autogenerated/w.cwl').read_text(encoding='utf-8'))['steps'][0]['run']

    for written in (tree, embedded):
        assert _entry(written) == SCRIPT
        assert written['outputs']['out']['type'] == 'File?'
        assert written['successCodes'] == [0, 1]


@pytest.mark.needs_cwltool
@pytest.mark.parametrize('runner', RUNNERS)
def test_a_packed_tool_and_one_as_written_run_wired_by_inference(tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
                                                                 runner: str) -> None:
    """`make`, packed, feeds `use`, as written, through the edge inference made; `use` prints `make`'s seed."""
    tools = tmp_path / 'tools'
    tools.mkdir()
    (tools / 'seed.txt').write_text('seed\n', encoding='utf-8')
    (tools / 'make.cwl').write_text(MAKE, encoding='utf-8')
    (tools / 'use.cwl').write_text(USE, encoding='utf-8')
    (tmp_path / 'w.wic').write_text(MAKE_USE, encoding='utf-8')
    monkeypatch.chdir(tmp_path)

    _run(tmp_path / 'w.wic', '--cwl_runner', runner)

    assert _produced() == ['seed\n', 'seed\n']


def _cat_config(directory: Path, config: str) -> Path:
    """A tool that prints its `config.txt` default, beside it in `directory`."""
    directory.mkdir(parents=True)
    (directory / 'config.txt').write_text(config, encoding='utf-8')
    tool = directory / 'echo.cwl'
    tool.write_text('''cwlVersion: v1.2
class: CommandLineTool
baseCommand: cat
stdout: out.txt
inputs:
  config: {type: File, default: {class: File, location: config.txt}, inputBinding: {position: 1}}
outputs:
  out: {type: stdout}
''', encoding='utf-8')
    return tool


@pytest.mark.needs_cwltool
@pytest.mark.parametrize('runner', RUNNERS)
def test_one_basepath_runs_two_echo_cwl_and_two_config_txt_with_provenance(tmp_path: Path,
                                                                           monkeypatch: pytest.MonkeyPatch,
                                                                           runner: str) -> None:
    """Two tools both called `echo.cwl`, each with its own `config.txt` default, run from one base path with
    provenance; each step prints its own tool's file."""
    _cat_config(tmp_path / 'dataA', 'config of A\n')
    _cat_config(tmp_path / 'dataB', 'config of B\n')
    (tmp_path / 'w.wic').write_text('steps:\n- id: a\n  run: dataA/echo.cwl\n- id: b\n  run: dataB/echo.cwl\n',
                                    encoding='utf-8')
    monkeypatch.chdir(tmp_path)

    _run(tmp_path / 'w.wic', '--cwl_runner', runner)

    assert _produced() == ['config of A\n', 'config of B\n']
    assert (Path('autogenerated') / 'provenance' / 'w').is_dir()


def _hashes(root: Path) -> dict[str, tuple[str, int]]:
    return {path.relative_to(root).as_posix(): (hashlib.sha256(path.read_bytes()).hexdigest(), path.stat().st_mtime_ns)
            for path in root.rglob('*') if path.is_file()}


@pytest.mark.needs_cwltool
def test_nothing_outside_the_basepath_is_written_or_changed(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Compiling and running packed tools reads the tools, scripts and data where they are and writes none of
    them, nor anything beside them."""
    project = tmp_path / 'project'
    _relcat(project / 'tools')
    _cat_config(project / 'dataA', 'config of A\n')
    (project / 'w.wic').write_text('steps:\n- id: relcat\n  run: tools/relcat.cwl\n- id: a\n  run: dataA/echo.cwl\n',
                                   encoding='utf-8')
    before = _hashes(project)
    (tmp_path / 'work').mkdir()
    monkeypatch.chdir(tmp_path / 'work')

    _run(project / 'w.wic')

    assert _hashes(project) == before
    assert _produced() == ['config of A\n', PRINTED]


@pytest.mark.fast
def test_packing_does_not_fetch_a_tools_schemas(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """`$schemas` stay as written: packing a tool never goes to the network for an ontology."""
    ontology = 'http://schemas.invalid/ontology.owl'
    fetched: list[str] = []
    fetch_text = DefaultFetcher.fetch_text

    def spy(self: DefaultFetcher, url: str, content_types: list[str] | None = None) -> str:
        fetched.append(url)
        return fetch_text(self, url, content_types)

    monkeypatch.setattr(DefaultFetcher, 'fetch_text', spy)
    tool = yaml.safe_load(PLAIN)
    tool['$schemas'] = [ontology]

    packed = utils_cwl.packed(yaml.safe_dump(tool), (tmp_path / 'echo.cwl').as_uri())

    assert ontology not in fetched
    assert packed['$schemas'] == [ontology]
