"""How a run is, and how it ended: what `Workflow.start()`'s handle reports.

The states and the status are those the compute service reports for a run (its contract 2.0.0),
so a run reads the same wherever it runs.
"""
from dataclasses import dataclass
from datetime import datetime
from enum import Enum
from pathlib import Path

from .wic_types import Json


class RunState(Enum):
    """Where a run is: the states the compute service reports, too.

    QUEUED and RUNNING may change; COMPLETED, FAILED and CANCELLED never do. A local
    run is never seen QUEUED: `Workflow.start()` returns once its runner has started,
    or could not be started.
    """

    QUEUED = 'QUEUED'
    RUNNING = 'RUNNING'
    COMPLETED = 'COMPLETED'
    FAILED = 'FAILED'
    CANCELLED = 'CANCELLED'


@dataclass(frozen=True, slots=True)
class RunResult:
    """How a run ended, and where it left its outputs and its log."""

    #: COMPLETED, FAILED or CANCELLED.
    state: RunState
    #: The workflow's output object, when the run COMPLETED; else None.
    outputs: Json | None
    #: The directory the runner moved the workflow's outputs to.
    outdir: Path
    #: Both of the runner's streams, as they came: `<basepath>/run_<workflow>.log`.
    log_file: Path
    #: The runner's exit code; 1 when it could not be started, negative when a signal ended it.
    exit_code: int


@dataclass(frozen=True, slots=True)
class StepStatus:
    """A step of the run that failed: how it ended, by exit status or by signal, when the runner said."""

    name: str
    state: RunState
    exit_status: int | None = None
    signal: str | None = None


@dataclass(frozen=True, slots=True)
class RunStatus:  # pylint: disable=too-many-instance-attributes
    """A run's status: where it is, and how it ended once it has."""

    #: Where the run is.
    state: RunState
    #: Why it is in its state, when there is something to say.
    reason: str | None = None
    #: When the runner started; None until then.
    started_at: datetime | None = None
    #: When the run ended; None until then.
    ended_at: datetime | None = None
    #: The steps that failed, in the order they failed.
    steps: tuple[StepStatus, ...] = ()
    #: The workflow's output object, once the run COMPLETED.
    outputs: Json | None = None
    #: The runner's exit code, once it exited.
    exit_code: int | None = None
