"""Running a workflow on this machine: the one path every front door takes.

A front door describes the run as a `backends.RunRequest`; `run` starts the runner its
backend launches as a child process, streams what it prints, and reports how it ended.
`start` returns the `Run` handle that `run` waits on.
"""
import enum
import os
import shlex
import signal
import stat
import subprocess as sub
import sys
import threading
import traceback
from contextlib import AbstractContextManager
from dataclasses import dataclass, field, replace
from datetime import datetime
from pathlib import Path
from collections.abc import Mapping
from typing import IO, Final

from sophios.wic_types import Json
from . import auto_gen_header, backends, preflight, realtime
from . import utils
from .backends import Backend, Launch, RunRequest
from .compute_request import ComputeRequest
from .plugins import logging_filters


@dataclass(frozen=True, slots=True)
class _CompiledWorkflowForCompute:
    name: str
    cwl_workflow: Json
    cwl_job_inputs: Json


def generate_run_script(cmdline: str) -> None:
    """Writes the command used to invoke the cwl-runner to run.sh
    Does not actually invoke ./run.sh

    Args:
        cmdline (str): The command line which invokes the cwl-runner
    """
    with open('run.sh', mode='w', encoding='utf-8') as f:
        f.write('#!/bin/bash -e\n')
        f.write(auto_gen_header)
        f.write(cmdline)
    # chmod +x run.sh
    st = os.stat('run.sh')
    os.chmod('run.sh', st.st_mode | stat.S_IEXEC)


class RunState(enum.Enum):
    """Where a run is. QUEUED and RUNNING may change; COMPLETED, FAILED and CANCELLED never do."""

    QUEUED = 'QUEUED'
    RUNNING = 'RUNNING'
    COMPLETED = 'COMPLETED'
    FAILED = 'FAILED'
    CANCELLED = 'CANCELLED'


@dataclass(frozen=True, slots=True)
class StepStatus:
    """A step of the run that failed: how it ended, by exit status or by signal, when the runner said."""

    name: str
    state: RunState
    exit_status: int | None = None
    signal: str | None = None


@dataclass(frozen=True, slots=True)
class RunStatus:  # pylint: disable=too-many-instance-attributes
    """A run's status.

    Attributes:
        state: Where the run is.
        reason: Why it is in its state, when there is something to say.
        started_at: When the runner started; None until then.
        ended_at: When the run ended; None until then.
        steps: The steps that failed, in the order they failed.
        outputs: The workflow's output object, once the run COMPLETED.
        exit_code: The runner's exit code, once it exited.
    """

    state: RunState
    reason: str | None = None
    started_at: datetime | None = None
    ended_at: datetime | None = None
    steps: tuple[StepStatus, ...] = ()
    outputs: Json | None = None
    exit_code: int | None = None


#: What an exit status means on its own (shell conventions), said after it.
_STATUS_MEANS: Final = {
    126: 'its command could not be executed (not executable, or built for another platform)',
    127: "its command was not found: the tool's baseCommand is not installed in its container image, "
         'or on this machine for a tool with none',
    137: 'it was killed, usually for running out of memory: give the container engine more memory '
         '(Docker Desktop: Settings, Resources), or the step less input',
}


def _failure_lines(failed: Mapping[str, int | str | None]) -> list[str]:
    """One line per failed step: who, how it ended, what that means, where its messages are."""
    if not failed:
        return ['Failure! Above, the first ERROR line, or the error the runner itself raised, says why the run failed.']
    lines = []
    for job, ended in failed.items():
        means: str | None
        if isinstance(ended, int):
            how, means = f'exited with status {ended}', _STATUS_MEANS.get(ended)
        elif ended is not None:
            how, means = f'was ended by {ended}', _STATUS_MEANS[137] if ended == 'SIGKILL' else None
        else:
            how, means = 'could not run, or its outputs could not be collected', None
        lines.append(f"Failure! step '{job}' {how}" + (f', which means {means}' if means else '')
                     + '. Its own messages are above, under its name.')
    return lines


def log_path(request: RunRequest) -> Path:
    """The file a run's two streams are written to, as they come: `<basepath>/run_<workflow>.log`."""
    return request.basepath / f'run_{request.workflow_name}.log'


#: How long a runner asked to stop has to stop the containers it started before it is killed.
STOP_GRACE_SECONDS: Final = 30.0


def _spawn(launch: Launch) -> sub.Popen[bytes]:
    """Start `launch` with its streams piped, in its own session, where a terminal's Ctrl-C does not reach it:
    Sophios stops it.

    The program is this Python environment's own when it has one (`preflight.runner_program`). POSIX
    only: no runner runs on Windows, and the pre-flight stops a run there first (`PLATFORM_LACKS_RUNNERS`).
    """
    if sys.platform == 'win32':
        raise AssertionError('the pre-flight stops a run on Windows')
    executable = preflight.runner_program(launch.argv[0])
    return sub.Popen(launch.argv, executable=executable, env=dict(launch.env), cwd=launch.cwd,
                     stdin=sub.DEVNULL, stdout=sub.PIPE, stderr=sub.PIPE, start_new_session=True)


@dataclass
class _Streams:
    """What the runner printed, as the run goes: echoed, written to the log, and its stdout kept."""

    log: IO[str]
    stdout: list[str] = field(default_factory=list)
    lock: threading.Lock = field(default_factory=threading.Lock)
    error: BaseException | None = None

    def pump(self, stream: IO[bytes], name: str) -> None:
        """Copy `stream` line by line to `sys.<name>` and to the log, until the runner closes it."""
        try:
            for raw in iter(stream.readline, b''):
                line = raw.decode('utf-8', errors='replace')
                echo = getattr(sys, name)
                echo.write(line)
                echo.flush()
                with self.lock:
                    self.log.write(line)
                    self.log.flush()
                    if name == 'stdout':
                        self.stdout.append(line)
        except Exception as error:  # pylint: disable=broad-exception-caught
            self.error = error
            for _ in iter(stream.readline, b''):   # keep draining, so the runner never blocks on a full pipe
                pass
        finally:
            stream.close()


class Run:  # pylint: disable=too-many-instance-attributes
    """A runner started as a child of this process, for as long as this process lives.

    `wait` blocks until it ends; `status`, `logs` and `outputs` say where it is; `cancel`
    stops it, as Ctrl-C during `wait` does.
    """

    def __init__(self, request: RunRequest, backend: Backend, launch: Launch,
                 plans: tuple[realtime.Plan, ...], analysis: realtime.Command | None) -> None:
        self.request = request
        self._backend = backend
        self._launch = launch
        self._plans = plans
        self._analysis = analysis
        self._status = RunStatus(RunState.QUEUED)
        self._process: sub.Popen[bytes] | None = None
        self._streams: _Streams | None = None
        self._readers: list[threading.Thread] = []
        self._watching: AbstractContextManager[None] | None = None
        self._settle_lock = threading.Lock()

    def start(self) -> None:
        """Start the runner, and the real-time analyses beside it.

        A runner that cannot be started is a FAILED run, its traceback in `error_<workflow>.txt`.
        """
        log = log_path(self.request)
        print('Running ' + shlex.join(self._launch.argv), flush=True)
        if self._analysis is not None:
            self._watching = realtime.watching(self._plans, Path(self.request.cachedir), self._analysis,
                                               env=self._launch.env)
            self._watching.__enter__()  # pylint: disable=unnecessary-dunder-call
        started_at = datetime.now()
        try:
            log.parent.mkdir(parents=True, exist_ok=True)
            self._streams = _Streams(log.open('w', encoding='utf-8'))
            self._process = _spawn(self._launch)
        except Exception as error:  # pylint: disable=broad-exception-caught
            self._unsupervised(error)
            if self._streams is not None:
                self._streams.log.close()
            self._stop_watching(cancelled=True)
            for line in _failure_lines({}):
                print(line, file=sys.stderr)
            self._status = RunStatus(RunState.FAILED, f'the runner could not be started: {error}', started_at,
                                     datetime.now(), exit_code=1)
            return
        self._status = RunStatus(RunState.RUNNING, started_at=started_at)
        try:
            assert self._process.stdout is not None and self._process.stderr is not None
            self._readers = [threading.Thread(target=self._streams.pump, args=(pipe, name), daemon=True)
                             for pipe, name in ((self._process.stdout, 'stdout'), (self._process.stderr, 'stderr'))]
            for reader in self._readers:
                reader.start()
        except KeyboardInterrupt:
            self.cancel()
            raise

    def _unsupervised(self, error: BaseException) -> None:
        """Say that Sophios could not start or follow the runner; keep the traceback in `error_<workflow>.txt`."""
        name = self.request.workflow_name
        print('Failed to execute', self.request.basepath / name)
        print(error)
        print(f'See error_{name}.txt for detailed technical information.')
        # Do not display a nasty stack trace to the user; hide it in a file.
        with open(f'error_{name}.txt', mode='w', encoding='utf-8') as f:
            traceback.print_exception(error, file=f)

    def _stop_watching(self, *, cancelled: bool) -> None:
        """End the real-time analyses: with their final runs, or at once when the run was cancelled."""
        watching, self._watching = self._watching, None
        if watching is None:
            return
        if cancelled:
            interrupt = KeyboardInterrupt()
            watching.__exit__(KeyboardInterrupt, interrupt, None)
        else:
            watching.__exit__(None, None, None)

    def _join_readers(self) -> None:
        for reader in self._readers:
            while reader.is_alive():
                reader.join(0.2)

    def wait(self) -> RunStatus:
        """Block until the run ends, and return its status.

        Raises:
            KeyboardInterrupt: On Ctrl-C, after the runner was stopped (`cancel`).
        """
        try:
            if self._process is not None:
                while self._process.poll() is None:
                    try:
                        self._process.wait(0.2)
                    except sub.TimeoutExpired:
                        pass
                self._join_readers()
        except KeyboardInterrupt:
            self.cancel()
            raise
        return self._settle()

    def status(self) -> RunStatus:
        """Where the run is now, without blocking. A runner that has exited, its output read to the end,
        is reported on, as `wait` does; while something it started still holds its output, it is RUNNING."""
        if self._process is not None and self._process.poll() is not None \
                and not any(reader.is_alive() for reader in self._readers):
            return self._settle()
        return self._status

    def logs(self) -> str:
        """Everything the runner has printed so far, both streams as they came."""
        try:
            return log_path(self.request).read_text(encoding='utf-8', errors='replace')
        except FileNotFoundError:
            return ''

    def outputs(self) -> Json | None:
        """The workflow's output object, once the run COMPLETED; else None."""
        return self.status().outputs

    def cancel(self) -> None:
        """Stop the runner and the analyses beside it, and mark the run CANCELLED.

        On POSIX the runner is sent SIGTERM, on which cwltool stops the containers it started,
        and killed with everything it started when it has not exited within `STOP_GRACE_SECONDS`.
        """
        # Held throughout, so a `wait` in another thread reports the run CANCELLED and not FAILED.
        with self._settle_lock:
            process = self._process
            if process is None or self._status.state not in (RunState.QUEUED, RunState.RUNNING):
                return
            if process.poll() is None:
                try:
                    process.send_signal(signal.SIGTERM)
                    process.wait(STOP_GRACE_SECONDS)
                except (sub.TimeoutExpired, KeyboardInterrupt):
                    self._kill(process)
                    process.wait()
            # Whatever the runner left in its session, which would also keep its streams open.
            self._kill(process)
            self._join_readers()
            if self._streams is not None:
                self._streams.log.close()
            self._stop_watching(cancelled=True)
            print("Interrupted; the runner's processes were terminated.", file=sys.stderr)
            self._status = RunStatus(RunState.CANCELLED, 'cancelled: the runner was stopped',
                                     self._status.started_at, datetime.now(), exit_code=process.returncode)

    @staticmethod
    def _kill(process: sub.Popen[bytes]) -> None:
        """Kill the runner and every process left in its session."""
        if sys.platform == 'win32':
            raise AssertionError('the pre-flight stops a run on Windows')
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except (ProcessLookupError, PermissionError):
            pass   # nothing is left in it

    def _settle(self) -> RunStatus:
        """Report on a runner that has exited, once: its outcome, the analyses' final runs, the summary lines."""
        with self._settle_lock:
            if self._status.state not in (RunState.QUEUED, RunState.RUNNING) or self._process is None:
                return self._status
            assert self._streams is not None
            exit_code = self._process.returncode
            self._streams.log.close()
            outcome = self._backend.outcome(self.request, exit_code, ''.join(self._streams.stdout),
                                            log_path(self.request))
            if outcome.summary is not None:
                print(f'Final output json metadata blob is in output_{self.request.workflow_name}.json')
            self._stop_watching(cancelled=False)
            if self._streams.error is not None:
                self._unsupervised(self._streams.error)
            steps = tuple(StepStatus(name, RunState.FAILED, ended if isinstance(ended, int) else None,
                                     ended if isinstance(ended, str) else None)
                          for name, ended in outcome.failed.items())
            if exit_code == 0:
                print(f'Success! Runner outputs are under {self._launch.outdir}/')
                self._status = RunStatus(RunState.COMPLETED, None, self._status.started_at, datetime.now(),
                                         steps, outcome.outputs, exit_code=0)
            else:
                for line in _failure_lines(outcome.failed):
                    print(line, file=sys.stderr)
                self._status = RunStatus(RunState.FAILED, f'the runner exited with status {exit_code}',
                                         self._status.started_at, datetime.now(), steps, exit_code=exit_code)
            return self._status


def _realtime(request: RunRequest, backend: Backend) -> tuple[RunRequest, realtime.Command | None]:
    """`request` ready to run its real-time analyses, and how each is run; or, when it cannot, a line saying
    why and no analyses. With analyses and no cachedir, cwltool caches in `cachedir/`."""
    if not request.realtime_plans:
        return request, None
    if request.generate_run_script:
        print('Real-time analysis runs only with --run_local; run.sh runs the workflow without it.',
              file=sys.stderr)
        return replace(request, realtime_plans=()), None
    command = backend.realtime_command(request)
    if isinstance(command, str):
        print(command, file=sys.stderr)
        return replace(request, realtime_plans=()), None
    if not request.cachedir:
        request = replace(request, cachedir='cachedir')
        print(f"Real-time analysis watches the runner's cache, so cwltool caches in {request.cachedir}/")
    return request, command


def start(request: RunRequest) -> Run:
    """Start the run `request` describes, and return its handle.

    Raises:
        ValueError: No backend runs the request's target with its runner, or it asks for a run script.
    """
    if request.generate_run_script:
        raise ValueError('a run script is written by run(), not started')
    backend = backends.backend(request.target, request.runner)
    request, analysis = _realtime(request, backend)
    handle = Run(request, backend, backend.launch(request), request.realtime_plans, analysis)
    handle.start()
    return handle


def run(request: RunRequest) -> int:
    """Run the workflow `request` describes and wait for it, or write run.sh for it.

    Returns:
        int: 0 on success, else the runner's exit code; 1 when the runner could not be started.

    Raises:
        KeyboardInterrupt: On Ctrl-C, after the runner was stopped.
        ValueError: No backend runs the request's target with its runner.
    """
    backend = backends.backend(request.target, request.runner)
    if request.generate_run_script:
        _realtime(request, backend)
        generate_run_script(shlex.join(backend.launch(request).argv))
        return 0  # Do not actually run
    status = start(request).wait()
    assert status.exit_code is not None
    return status.exit_code


def run_compute(workflow_name: str, workflow: Json, workflow_inputs: Json,
                submit_url: str) -> int | None:
    """Submit a compiled workflow to compute-slurm.

    Args:
        workflow_name (str): The name of the workflow.
        workflow (Json): The compiled CWL workflow.
        workflow_inputs (Json): The inputs for compiled CWL workflow.
        submit_url (str): URL of Compute where the job is to be submitted.
    Returns:
        int | None: The return value indicating if submission succeeded (`0`) or not.
    """
    now = datetime.now()
    date_time = now.strftime("%Y_%m_%d_%H.%M.%S")
    jobid = workflow_name + '__' + str(date_time) + '__'
    compute_request = ComputeRequest(
        _CompiledWorkflowForCompute(workflow_name, workflow, workflow_inputs),
        workflow_id=jobid,
    )

    # sanity check if the string has the form of an URL
    if not utils.is_valid_url(submit_url):
        print("Ill-formed URL string detected! Please provide a valid URL")
        return 1

    submission = compute_request.submit(
        submit_url,
        log_path=Path(f'compute_logs_{jobid}.txt'),
    )
    return submission.exit_code


def cwltool_main() -> int:
    """Another entrypoint for filtering regular logged output.

    `cwltool.main.run`, not `main`: it stops the containers the run started when the run is sent SIGTERM.
    """
    import cwltool.main  # pylint: disable=import-outside-toplevel
    logging_filters()
    retval: int = cwltool.main.run(sys.argv[1:])
    return retval


def cwltool_main_pf() -> int:
    """Another entrypoint for filtering partial failures logged output"""
    import cwltool.main  # pylint: disable=import-outside-toplevel
    logging_filters(True)
    retval: int = cwltool.main.run(sys.argv[1:])
    return retval
