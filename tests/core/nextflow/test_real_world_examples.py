"""Real-world pipeline shapes written in Sophios and run by Nextflow.

Each test builds one shape in Python, writes the generated pipeline, runs it
under the pinned Nextflow, and compares the one report it writes with a literal.
The shapes: a registry of sources that decides the tasks; one step per provider
over (year, month) pairs; per-year downloads feeding aggregates; checks behind
every source merged into one input; boolean flags set in the params file and on
the command line; and a preflight step whose output reaches every scattered task.

Every task list is enumerated in Python when the workflow is built, and every
dependency is data that a task reads.
"""

# pylint: disable=missing-function-docstring

import json
from pathlib import Path
from typing import Any

import pytest

from sophios.api.python.tool_builder import CommandLineTool, Input, Inputs, Output, Outputs, cwl
from sophios.api.python.workflow import Step, StepInput, Workflow
from sophios.nf_types import NfWorkflowInputConnection

from .testkit import execute_nextflow

# One line per region the task covers.
_EMIT = 'source=$1; year=$2; shift 2; for region in "$@"; do printf "%s %s %s\\n" "$source" "$year" "$region"; done'
# One line, from the download among the arguments that is this task's year.
_FROM_DOWNLOAD = ('year=$1; region=$2; shift 2; for data in "$@"; do read -r kind of < "$data"; '
                  '[ "$of" = "$year" ] && printf "alpha %s %s %s\\n" "$year" "$region" "$kind"; done; true')
# A source's line count, then its lines, across every shard of every source.
_CHECK = ('source=$1; shift; lines=$(cat "$@" </dev/null | grep "^$source " || true); '
          'printf "%s %s\\n" "$source" "$(printf "%s" "$lines" | grep -c . || true)"; '
          '[ -z "$lines" ] || printf "%s\\n" "$lines"')


def _sh(name: str, script: str, output: str, **inputs: Input) -> CommandLineTool:
    """A tool that runs `script` under sh with its inputs as positional parameters, and keeps stdout."""
    outputs = Outputs(**{output: Output(cwl.file, glob=f"{output}.txt")})
    tool = CommandLineTool(name, Inputs(**inputs), outputs)
    return tool.base_command("sh", "-c", script, name).stdout(f"{output}.txt")


def _report(**files: Any) -> Step:
    """One report: every file of each input, in input order."""
    tool = CommandLineTool(
        "report", Inputs(**{name: Input(cwl.array(cwl.file), position=index) for index, name in enumerate(files, 1)}),
        Outputs(report=Output(cwl.file, glob="report.txt")),
    ).base_command("cat").stdout("report.txt")
    report = Step(tool, step_name="report")
    for name, source in files.items():
        setattr(report.inputs, name, source)
    return report


def _run_report(workflow: Workflow, directory: Path, *arguments: str) -> str:
    """Write ``workflow``'s pipeline, run it with ``arguments``, and read the one report it writes."""
    workflow.to_nextflow(directory)
    run = execute_nextflow(directory, *arguments)
    assert run.returncode == 0, f"stdout:\n{run.stdout}\nstderr:\n{run.stderr}"
    reports = list((directory / "work").rglob("report.txt"))
    assert len(reports) == 1, f"no single report was written:\n{run.stdout}"
    return reports[0].read_text(encoding="utf-8")


def _regions_for(year: int) -> list[str]:
    """Zone regions change version with the year; grid regions do not."""
    return [f"zone_v{1 if year <= 2001 else 2}", "grid_5"]


# Modelled on Matthew's repo, main.nf.
@pytest.mark.nextflow
@pytest.mark.serial
def test_a_registry_decides_each_sources_tasks(tmp_path: Path) -> None:
    """(year, region) pairs; a crossproduct over a fixed list; one task per year; a source with no tasks."""
    per_region = _sh("per_region", _EMIT, "shard", source=Input(cwl.string, position=1),
                     year=Input(cwl.int, position=2), region=Input(cwl.string, position=3))

    download = Step(_sh("download", 'printf "data %s\\n" "$1"', "data", year=Input(cwl.int, position=1)),
                    step_name="download")
    download.inputs.year = [2001, 2002]
    download.scatter_on(download.inputs.year)

    # alpha: (year, region) pairs, the regions depending on the year; each task is given
    # every download and keeps the one for its own year.
    pairs = [(year, region) for year in (2001, 2002) for region in _regions_for(year)]
    alpha = Step(_sh("from_download", _FROM_DOWNLOAD, "shard", year=Input(cwl.int, position=1),
                     region=Input(cwl.string, position=2), downloads=Input(cwl.array(cwl.file), position=3)),
                 step_name="alpha")
    alpha.inputs.year = [year for year, _ in pairs]
    alpha.inputs.region = [region for _, region in pairs]
    alpha.inputs.downloads = download.outputs.data
    alpha.scatter_on(alpha.inputs.year, alpha.inputs.region)

    # beta: a fixed region list crossed with the years.
    beta = Step(per_region, step_name="beta")
    beta.inputs.source = "beta"
    beta.inputs.year = [2002, 2003]
    beta.inputs.region = ["grid_3", "grid_5"]
    beta.scatter_on(beta.inputs.year, beta.inputs.region, method="flat_crossproduct")

    # delta: one task per year over all of that year's regions.
    delta = Step(_sh("per_year", _EMIT, "shard", source=Input(cwl.string, position=1),
                     year=Input(cwl.int, position=2), regions=Input(cwl.array(cwl.string), position=3)),
                 step_name="delta")
    delta.inputs.source = "delta"
    delta.inputs.year = [2002, 2003]
    delta.inputs.regions = _regions_for(2002)
    delta.scatter_on(delta.inputs.year)

    # gamma: no years, so no tasks; its check still runs.
    gamma = Step(per_region, step_name="gamma")
    gamma.inputs.source = "gamma"
    gamma.inputs.year = []
    gamma.inputs.region = []
    gamma.scatter_on(gamma.inputs.year, gamma.inputs.region)

    sources = {"alpha": alpha, "beta": beta, "delta": delta, "gamma": gamma}
    check = Step(_sh("check", _CHECK, "check", source=Input(cwl.string, position=1),
                     **{f"shards_{name}": Input(cwl.array(cwl.file), position=2) for name in sources}),
                 step_name="check")
    check.inputs.source = list(sources)
    for name, source in sources.items():
        setattr(check.inputs, f"shards_{name}", source.outputs.shard)
    check.scatter_on(check.inputs.source)

    workflow = Workflow([download, *sources.values(), check, _report(checks=check.outputs.check)], "registry")
    assert _run_report(workflow, tmp_path) == """\
alpha 4
alpha 2001 zone_v1 data
alpha 2001 grid_5 data
alpha 2002 zone_v2 data
alpha 2002 grid_5 data
beta 4
beta 2002 grid_3
beta 2002 grid_5
beta 2003 grid_3
beta 2003 grid_5
delta 4
delta 2002 zone_v2
delta 2002 grid_5
delta 2003 zone_v2
delta 2003 grid_5
gamma 0
"""


def _months(year: int, month: int, count: int) -> list[tuple[int, int]]:
    """``count`` consecutive (year, month) pairs from ``year``-``month``."""
    return [(year + (month - 1 + i) // 12, (month - 1 + i) % 12 + 1) for i in range(count)]


# Modelled on Matthew's repo, satellite.nf.
@pytest.mark.nextflow
@pytest.mark.serial
def test_one_step_per_provider_runs_one_task_per_month(tmp_path: Path) -> None:
    """Each provider's step is a dotproduct over (year, month) pairs enumerated in Python."""
    per_month = _sh("per_month", 'printf "%s %04d-%02d\\n" "$1" "$2" "$3"', "shard",
                    provider=Input(cwl.string, position=1), year=Input(cwl.int, position=2),
                    month=Input(cwl.int, position=3))
    steps = {}
    for provider, pairs in {"provider_a": _months(2001, 11, 4), "provider_b": _months(2002, 6, 2)}.items():
        step = steps[provider] = Step(per_month, step_name=provider)
        step.inputs.provider = provider
        step.inputs.year = [year for year, _ in pairs]
        step.inputs.month = [month for _, month in pairs]
        step.scatter_on(step.inputs.year, step.inputs.month)
    report = _report(**{name: step.outputs.shard for name, step in steps.items()})
    assert _run_report(Workflow([*steps.values(), report], "months"), tmp_path) == """\
provider_a 2001-11
provider_a 2001-12
provider_a 2002-01
provider_a 2002-02
provider_b 2002-06
provider_b 2002-07
"""


# Modelled on Matthew's repo, nlcd.nf.
@pytest.mark.nextflow
@pytest.mark.serial
def test_per_year_downloads_feed_scattered_aggregates(tmp_path: Path) -> None:
    """A year's download reaches every task of its year's aggregate; another aggregate waits on every download."""
    download_tool = _sh("download", 'printf "data %s\\n" "$1"', "data", year=Input(cwl.int, position=1))
    aggregate_tool = _sh("aggregate", 'read -r kind of < "$2"; printf "%s %s %s\\n" "$of" "$1" "$kind"', "summary",
                         region=Input(cwl.string, position=1), data=Input(cwl.file, position=2))
    downloads, aggregates = [], {}
    for year in (2001, 2002):
        download = Step(download_tool, step_name=f"download_{year}")
        download.inputs.year = year
        aggregate = aggregates[f"aggregate_{year}"] = Step(aggregate_tool, step_name=f"aggregate_{year}")
        aggregate.inputs.region = _regions_for(year)
        aggregate.inputs.data = download.outputs.data
        aggregate.scatter_on(aggregate.inputs.region)
        downloads.append(download)

    # Every task is given every download and keeps the one for its own year.
    pairs = [(year, region) for year in (2001, 2002) for region in _regions_for(year)]
    across = Step(_sh("from_download", _FROM_DOWNLOAD, "summary", year=Input(cwl.int, position=1),
                      region=Input(cwl.string, position=2), downloads=Input(cwl.array(cwl.file), position=3)),
                  step_name="across")
    across.inputs.year = [year for year, _ in pairs]
    across.inputs.region = [region for _, region in pairs]
    across.inputs.downloads = StepInput(source=[download.outputs.data for download in downloads],
                                        link_merge="merge_flattened")
    across.scatter_on(across.inputs.year, across.inputs.region)

    report = _report(**{name: step.outputs.summary for name, step in {**aggregates, "across": across}.items()})
    assert _run_report(Workflow([*downloads, *aggregates.values(), across, report], "aggregates"), tmp_path) == """\
2001 zone_v1 data
2001 grid_5 data
2002 zone_v2 data
2002 grid_5 data
alpha 2001 zone_v1 data
alpha 2001 grid_5 data
alpha 2002 zone_v2 data
alpha 2002 grid_5 data
"""


# Modelled on Matthew's repo, quality.nf.
@pytest.mark.nextflow
@pytest.mark.serial
def test_a_check_per_source_reads_every_sources_shards_merged_into_one_input(tmp_path: Path) -> None:
    """Scattered sources, one with no tasks, merge_flattened into one input of a scattered check."""
    per_region = _sh("per_region", _EMIT, "shard", source=Input(cwl.string, position=1),
                     year=Input(cwl.int, position=2), region=Input(cwl.string, position=3))
    sources = {name: Step(per_region, step_name=name) for name in ("alpha", "beta", "gamma")}
    for name, source in sources.items():
        source.inputs.source = name
    alpha, beta, gamma = sources.values()
    alpha.inputs.year, alpha.inputs.region = [2001, 2002], ["zone_v1", "zone_v2"]
    beta.inputs.year, beta.inputs.region = [2002, 2003], ["grid_3", "grid_5"]
    gamma.inputs.year, gamma.inputs.region = [], []
    alpha.scatter_on(alpha.inputs.year, alpha.inputs.region)
    beta.scatter_on(beta.inputs.year, beta.inputs.region, method="flat_crossproduct")
    gamma.scatter_on(gamma.inputs.year, gamma.inputs.region)

    check = Step(_sh("check", _CHECK, "check", source=Input(cwl.string, position=1),
                     shards=Input(cwl.array(cwl.file), position=2)), step_name="check")
    check.inputs.source = list(sources)
    check.inputs.shards = StepInput(source=[source.outputs.shard for source in sources.values()],
                                    link_merge="merge_flattened")
    check.scatter_on(check.inputs.source)

    workflow = Workflow([*sources.values(), check, _report(checks=check.outputs.check)], "checks")
    assert _run_report(workflow, tmp_path) == """\
alpha 2
alpha 2001 zone_v1
alpha 2002 zone_v2
beta 4
beta 2002 grid_3
beta 2002 grid_5
beta 2003 grid_3
beta 2003 grid_5
gamma 0
"""


# Modelled on Matthew's repo, quality.nf.
@pytest.mark.nextflow
@pytest.mark.serial
def test_a_command_line_false_overrides_a_flag_the_params_file_sets(tmp_path: Path) -> None:
    """The params file sets strict and verbose and clears fast; the command line clears verbose."""
    flags = Step(_sh("flags", 'printf "%s\\n" "$@"', "report", strict=Input(cwl.boolean, position=1, flag="--strict"),
                     verbose=Input(cwl.boolean, position=2, flag="--verbose"),
                     fast=Input(cwl.boolean, position=3, flag="--fast")), step_name="flags")
    flags.inputs.strict, flags.inputs.verbose, flags.inputs.fast = True, True, False
    workflow = Workflow([flags], "flags")
    verbose = next(connection.from_port for connection in workflow.compile(target="nextflow").connections
                   if isinstance(connection, NfWorkflowInputConnection) and connection.to_port == "verbose")
    assert _run_report(workflow, tmp_path, f"--{verbose}", "false") == "--strict\n"
    assert json.loads((tmp_path / "nextflow_params.json").read_text(encoding="utf-8"))[verbose] is True


# Modelled on Matthew's repo, main.nf.
@pytest.mark.nextflow
@pytest.mark.serial
def test_a_preflight_steps_output_reaches_every_scattered_task(tmp_path: Path) -> None:
    preflight = Step(_sh("preflight", 'printf "checked %s\\n" "$1"', "status", target=Input(cwl.string, position=1)),
                     step_name="preflight")
    preflight.inputs.target = "grid_5"
    work = Step(_sh("work", 'read -r state < "$2"; printf "%s %s\\n" "$1" "$state"', "result",
                    tile=Input(cwl.string, position=1), status=Input(cwl.file, position=2)), step_name="work")
    work.inputs.tile = ["tile_1", "tile_2", "tile_3"]
    work.inputs.status = preflight.outputs.status
    work.scatter_on(work.inputs.tile)
    workflow = Workflow([preflight, work, _report(results=work.outputs.result)], "preflight")
    assert _run_report(workflow, tmp_path) == """\
tile_1 checked grid_5
tile_2 checked grid_5
tile_3 checked grid_5
"""
