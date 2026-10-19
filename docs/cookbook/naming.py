"""Naming conventions: the report should read the summary, not the details."""
from pathlib import Path

from sophios.api.python.workflow import Step, Workflow


#: The tools that ship with Sophios, where `sophios --generate_config` put them.
ADAPTERS = Path.home() / 'wic' / 'cwl_adapters'


def workflow() -> Workflow:
    """Write a text, summarize it, and report on it, every connection left to inference."""
    echo = Step(clt_path=ADAPTERS / 'echo.cwl')
    echo.inputs.message = 'The quick brown fox jumps over the lazy dog'
    summarize = Step(clt_path=Path(__file__).with_name('summarize.cwl'))
    report = Step(clt_path=Path(__file__).with_name('report.cwl'))
    return Workflow([echo, summarize, report], 'naming')


if __name__ == '__main__':
    workflow().run()
