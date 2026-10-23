"""What each backend gives a run: the runner's command line, and what it reads back from the runner."""
from pathlib import Path
from typing import Any

import pytest

from sophios import backends
from sophios.backends import RunRequest


def _argv(request: RunRequest) -> list[str]:
    return list(backends.backend(request.target, request.runner).launch(request).argv)


@pytest.mark.fast
@pytest.mark.parametrize('runner', ['cwltool', 'toil-cwl-runner'])
def test_a_given_outdir_is_used_once_and_clean_follows_the_runner(tmp_path: Path, runner: str) -> None:
    outdir = tmp_path / 'workflow_output'
    cmd = _argv(RunRequest('wf', tmp_path / 'exec', runner=runner, outdir=str(outdir)))
    assert cmd.count('--outdir') == 1
    assert cmd[cmd.index('--outdir') + 1] == str(outdir.resolve())
    if runner == 'toil-cwl-runner':
        assert cmd[cmd.index('--clean') + 1] == 'always'
    else:
        assert '--clean' not in cmd


@pytest.mark.fast
@pytest.mark.parametrize(('runner', 'name'), [('cwltool', 'cwltool'), ('toil-cwl-runner', 'toil')])
def test_two_workflows_started_in_the_same_second_get_their_own_outdirs(
        tmp_path: Path, monkeypatch: pytest.MonkeyPatch, runner: str, name: str) -> None:
    """The default output directory is named after the workflow, as provenance, the summary and the
    job store are, so two workflows started in the same second from one basepath never share it."""
    class Frozen(backends.datetime):
        @classmethod
        def now(cls, tz: Any = None) -> 'Frozen':
            return cls(2026, 10, 9, 12, 2, 16)

    monkeypatch.setattr(backends, 'datetime', Frozen)
    base = tmp_path.resolve()
    outdirs = []
    for workflow in ('bbbc', 'python_script_demo'):
        cmd = _argv(RunRequest(workflow, base, runner=runner))
        outdirs.append(cmd[cmd.index('--outdir') + 1])
    assert outdirs == [f'{base}/outdir_{name}_bbbc_2026_10_09_12.02.16',
                       f'{base}/outdir_{name}_python_script_demo_2026_10_09_12.02.16']


@pytest.mark.fast
def test_cwltool_is_quiet_unless_told_otherwise_and_toil_never_is(tmp_path: Path) -> None:
    """`--quiet` is on for callers that say nothing and off when `quiet=False`, so `--debug` can be heard."""
    assert '--quiet' in _argv(RunRequest('wf', tmp_path))
    assert '--quiet' not in _argv(RunRequest('wf', tmp_path, quiet=False))
    assert '--quiet' not in _argv(RunRequest('wf', tmp_path, runner='toil-cwl-runner'))


@pytest.mark.fast
@pytest.mark.parametrize(('cachedir', 'given'), [('', None), ('cache', 'cache')])
def test_cwltool_gets_the_cachedir_it_was_given_made_absolute(monkeypatch: pytest.MonkeyPatch, tmp_path: Path,
                                                              cachedir: str, given: str | None) -> None:
    monkeypatch.chdir(tmp_path)
    cmd = _argv(RunRequest('wf', tmp_path, cachedir=cachedir))
    if given is None:
        assert '--cachedir' not in cmd
    else:
        assert cmd[cmd.index('--cachedir') + 1] == str(tmp_path / given)


@pytest.mark.fast
@pytest.mark.parametrize(('passthrough', 'provenance'), [
    ((), True), (('--batchSystem', 'slurm'), False), (('--batchSystem=slurm',), True)])
def test_toil_keeps_no_provenance_when_an_argument_is_slurm(tmp_path: Path, passthrough: tuple[str, ...],
                                                            provenance: bool) -> None:
    """Toil writes no provenance under Slurm: an argument that is exactly `slurm` drops the flag."""
    cmd = _argv(RunRequest('wf', tmp_path, runner='toil-cwl-runner', passthrough_args=passthrough))
    assert ('--provenance' in cmd) is provenance


@pytest.mark.fast
def test_the_runner_gets_this_processs_environment_with_the_users_on_top(
        monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.setenv('SOPHIOS_TEST_KEPT', 'kept')
    monkeypatch.setenv('SOPHIOS_TEST_REPLACED', 'old')
    launch = backends.backend('cwl', 'cwltool').launch(
        RunRequest('wf', tmp_path, env={'SOPHIOS_TEST_REPLACED': 'new'}))
    assert launch.env['SOPHIOS_TEST_KEPT'] == 'kept' and launch.env['SOPHIOS_TEST_REPLACED'] == 'new'


@pytest.mark.fast
def test_an_unknown_runner_is_refused_naming_the_runners() -> None:
    with pytest.raises(ValueError, match="No backend runs cwl with 'cwl-runner'; the runners are cwltool"):
        backends.backend('cwl', 'cwl-runner')


@pytest.mark.fast
def test_failed_steps_are_read_from_the_runners_log_lines() -> None:
    lines = ['WARNING [job a] exited with status: 3', 'WARNING [job a] completed permanentFail',
             'INFO [job b] completed success',
             '\x1b[1;30mWARNING\x1b[0m \x1b[33m[job c] was terminated by signal: SIGKILL\x1b[0m',
             '\x1b[1;30mWARNING\x1b[0m \x1b[33m[job c] completed permanentFail\x1b[0m',
             '[2026-10-22 10:00:00] WARNING [job d] completed permanentFail']
    assert backends.failed_steps(lines) == {'a': 3, 'c': 'SIGKILL', 'd': None}


@pytest.mark.fast
def test_toil_writes_no_summary_when_it_printed_no_output_object(tmp_path: Path) -> None:
    request = RunRequest('wf', tmp_path, runner='toil-cwl-runner')
    backend = backends.backend('cwl', 'toil-cwl-runner')
    outcome = backend.outcome(request, 0, 'no object here\n', tmp_path / 'absent.log')
    assert outcome.outputs is None and outcome.summary is None
    assert not (tmp_path / 'output_wf.json').exists()
