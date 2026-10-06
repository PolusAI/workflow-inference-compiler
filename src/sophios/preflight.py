"""What a compile or a run needs from this machine, checked before the work starts.

The pre-flight is for local runs only: `sophios --run_local`, `--generate_run_script`
and `--check`, a plain CWL workflow run the same way, and `Workflow.run()`. It
asks about this machine's container engine, paths and programs, so a compute
submission, which runs elsewhere, never goes through it; what a submission checks
is its payload (the document compiles, the inputs are well formed).

Each check returns diagnostics rather than raising, so a caller reports every
problem at once, each on one line: what is wrong, then what to do.
"""
import os
import subprocess as sub
import sys
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Final

import cwl_utils.parser as cwl
from cwl_utils.docker_extract import traverse

from .lang.diagnostics import Diagnostic, Severity, SophiosError
from .lang.error_codes import SophiosErrorCode

#: Engines `<engine> info` reaches: it needs the engine's daemon or machine, and no container or network.
DOCKER_LIKE: Final = ('docker', 'podman')
#: More Docker Desktop processes than this can make a run hang.
MAX_DOCKER_PROCESSES: Final = 1000
#: Where the docker CLI looks for its daemon when `DOCKER_HOST` does not say.
DEFAULT_DOCKER_SOCKET: Final = Path('/var/run/docker.sock')


@dataclass(frozen=True, slots=True)
class RunSettings:
    """What a run was asked to do, as far as this machine is concerned."""

    container_engine: str
    pull_dir: str
    ignore_install: bool = False
    ignore_processes: bool = False
    #: Each directory the run writes into: the path, what it holds, what to do instead (see `unwritable`).
    writes: tuple[tuple[Path, str, str], ...] = ()


@dataclass(frozen=True, slots=True)
class Needs:
    """What running some CWL documents needs from this machine."""

    #: The documents the runner gets: the workflow first, then any real-time analyses.
    documents: tuple[Path, ...]
    #: Every DockerRequirement, as a requirement or a hint, of every process the documents run.
    containers: tuple[cwl.DockerRequirement, ...]


def needs(documents: Sequence[Path]) -> Needs:
    """Read what running `documents` needs, as cwl-docker-extract reads it (no `$schemas` are fetched)."""
    containers = [requirement for document in documents
                  for requirement in traverse(cwl.load_document_by_uri(str(document)))]
    return Needs(tuple(documents), tuple(containers))


def check(found: Needs, settings: RunSettings) -> None:
    """Raise every problem a run would hit on this machine, before anything is pulled.

    Raises:
        SophiosError: One diagnostic per problem.
    """
    problems = [problem for directory, holds, remedy in settings.writes
                if (problem := unwritable(directory, holds, remedy)) is not None]
    if found.containers:
        problems += _engine_problems(found, settings)
    if problems:
        raise SophiosError(problems)


def pull(found: Needs, settings: RunSettings) -> None:
    """Make every image the run uses available."""
    if not found.containers:
        return
    # cwl-docker-extract recursively `docker pull`s all images in all subworkflows.
    # This is important because cwltool only uses `docker run` when executing
    # workflows, and if there is a local image available,
    # `docker run` will NOT query the remote repository for the latest image!
    # cwltool has a --force-docker-pull option, but this may cause multiple pulls in parallel.
    for document in found.documents:
        if settings.container_engine == 'singularity':
            cmd = ['cwl-docker-extract', '-s', '--dir', settings.pull_dir, str(document)]
        else:
            cmd = ['cwl-docker-extract', '--force-download', str(document)]
        sub.run(cmd, check=True)


def prepare(documents: Sequence[Path], settings: RunSettings) -> None:
    """Check what running `documents` needs from this machine, then pull their images.

    The one call every local run makes before it starts: the CLI's compiled and
    plain-CWL runs and `Workflow.run()`.

    Raises:
        SophiosError: Every problem found, before anything is pulled.
    """
    found = needs(documents)
    check(found, settings)
    pull(found, settings)


def unwritable(directory: Path, holds: str, remedy: str) -> Diagnostic | None:
    """A `wic021` error when this user cannot create `directory` or write into it, else None.

    The nearest part of the path that exists decides: a directory that does not
    exist yet is fine when it can be created.

    Args:
        directory (Path): Where Sophios is about to write.
        holds (str): What it writes there, for the message: 'the compiled workflow'.
        remedy (str): What to do instead: 'run Sophios from a directory you can write to'.

    Returns:
        Diagnostic | None: The error, or None when `directory` can be written.
    """
    target = directory.absolute()
    existing = next(path for path in (target, *target.parents) if path.exists())
    if existing.is_dir() and os.access(existing, os.W_OK | os.X_OK):
        return None
    why = 'is not writable by you' if existing.is_dir() else 'is a file'
    return Diagnostic(Severity.ERROR, SophiosErrorCode.DIRECTORY_NOT_WRITABLE,
                      f'Sophios writes {holds} to {target}, but {existing} {why}: {remedy}.')


def _engine_problems(found: Needs, settings: RunSettings) -> list[Diagnostic]:
    """What stops the engine running the containers `found` needs; one diagnostic each.

    Only signals that do not depend on wording are classified: whether the engine
    program exists, the exit status of `<engine> info`, and whether the socket it
    would use exists and is accessible. What the engine said is quoted, never read.
    """
    engine = settings.container_engine
    probe = [engine, 'info'] if engine in DOCKER_LIKE else [engine, '--version']
    try:
        proc = sub.run(probe, check=False, stdout=sub.DEVNULL, stderr=sub.PIPE)
    except FileNotFoundError:
        if settings.ignore_install:
            return []
        return [_engine_error(f'{engine} is not installed (it is not on PATH), and this workflow runs tools in '
                              f'containers ({_images(found)}): install {engine}, or name an installed engine '
                              'with --container_engine')]
    problems = []
    if proc.returncode != 0 and not settings.ignore_install:
        problems.append(_unusable_engine(engine, proc.returncode, proc.stderr.decode('utf-8', errors='replace')))
    if engine in DOCKER_LIKE and sys.platform != 'win32' and not settings.ignore_processes:
        count = _docker_processes()
        if count > MAX_DOCKER_PROCESSES:
            problems.append(_engine_error(
                f'{count} docker processes are running, and more than {MAX_DOCKER_PROCESSES} can make a run '
                'hang: quit Docker and start it again (sudo pkill com.docker && sudo pkill Docker), or pass '
                '--ignore_docker_processes'))
    return problems


def _unusable_engine(engine: str, status: int, said: str) -> Diagnostic:
    """The engine is installed but its probe exited non-zero: its socket is not there, or not yours, or it said why.

    A socket is called missing only when `DOCKER_HOST` names it: without that, the CLI may reach its
    daemon through a context, and the default path says nothing.
    """
    line = next((text.strip() for text in said.splitlines() if text.strip()), 'it said nothing')
    socket = _docker_socket() if engine == 'docker' else None
    if socket is not None and 'DOCKER_HOST' in os.environ and not socket.exists():
        return _engine_error(f'{engine} is installed, but its engine is not reachable: the socket {socket} does '
                             f'not exist ({line}): start the engine, then run again')
    if socket is not None and not os.access(socket, os.R_OK | os.W_OK):
        return _engine_error(f'{engine} is installed, but you may not use its socket {socket} ({line}): add your '
                             'user to the docker group (sudo usermod -aG docker $USER), then log out and back in')
    return _engine_error(f'{engine} is installed, but `{engine} info` exited with status {status} ({line}): '
                         f'run `{engine} info` to see why, and start the engine if it is stopped')


def _docker_socket() -> Path | None:
    """The unix socket the docker CLI connects to: `DOCKER_HOST`'s, else the default; None for a remote host."""
    host = os.environ.get('DOCKER_HOST')
    if host is None:
        return DEFAULT_DOCKER_SOCKET
    return Path(host.removeprefix('unix://')) if host.startswith('unix://') else None


def _images(found: Needs) -> str:
    """The first image the run uses, and how many others: `docker.io/bash:4.4 and 2 more`."""
    names = list(dict.fromkeys(req.dockerPull or req.dockerImageId or 'an image built by the run'
                               for req in found.containers))
    return names[0] if len(names) == 1 else f'{names[0]} and {len(names) - 1} more'


def _docker_processes() -> int:
    """How many Docker Desktop processes run (it names them com.docker.*); 0 without pgrep."""
    try:
        listed = sub.run(['pgrep', 'com.docker'], check=False, stdout=sub.PIPE, stderr=sub.DEVNULL)
    except FileNotFoundError:
        return 0
    return len(listed.stdout.split())


def _engine_error(message: str) -> Diagnostic:
    return Diagnostic(Severity.ERROR, SophiosErrorCode.CONTAINER_ENGINE_UNAVAILABLE, f'{message}.')
