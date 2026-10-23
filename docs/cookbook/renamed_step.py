"""Renamed steps: one tool used twice, each step under its own name."""
from pathlib import Path

from sophios.api.python.workflow import Step, Workflow


#: The tools that ship with Sophios, where `sophios --generate_config` put them.
ADAPTERS = Path.home() / 'wic' / 'cwl_adapters'


def workflow() -> Workflow:
    """Run ``echo.cwl`` twice, as the steps ``hello`` and ``goodbye``."""
    hello = Step(clt_path=ADAPTERS / 'echo.cwl', step_name='hello')
    hello.inputs.message = 'Hello'
    goodbye = Step(clt_path=ADAPTERS / 'echo.cwl', step_name='goodbye')
    goodbye.inputs.message = 'Goodbye'
    return Workflow([hello, goodbye], 'renamed_step')


if __name__ == '__main__':
    workflow().run()
