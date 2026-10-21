"""Channels and pipeline safety under installed, pinned Nextflow.

design_docs/sophios_nextflow_backend.md §6, Inputs and channels and Topology;
§7 for the stub every process carries.
"""

# pylint: disable=missing-function-docstring

from pathlib import Path

import pytest

from sophios.input_output_nf import NF_STUB_MESSAGE, write_nextflow_artifacts

from .testkit import execute_nextflow, runtime_workflow


@pytest.mark.nextflow
@pytest.mark.serial
def test_a_dry_run_stops_at_its_first_task_and_writes_nothing(tmp_path: Path) -> None:
    write_nextflow_artifacts(runtime_workflow(), tmp_path)
    result = execute_nextflow(tmp_path, "-stub-run")
    assert result.returncode == 1, result.stdout + result.stderr
    assert NF_STUB_MESSAGE in result.stdout + result.stderr
    assert not list((tmp_path / "work").rglob("message.txt"))
    assert not list((tmp_path / "work").rglob("copy.txt"))
