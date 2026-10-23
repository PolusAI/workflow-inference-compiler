"""A named workflow output."""
from pathlib import Path

from sophios.api.python.workflow import Step, Workflow


#: The tools that ship with Sophios, where `sophios --generate_config` put them.
ADAPTERS = Path.home() / 'wic' / 'cwl_adapters'


def workflow() -> Workflow:
    """Echo a greeting and name the result ``greeting``."""
    echo = Step(clt_path=ADAPTERS / 'echo.cwl')
    echo.inputs.message = 'Hello'
    flow = Workflow([echo], 'named_output')
    flow.outputs.greeting = echo.outputs.stdout
    return flow


if __name__ == '__main__':
    workflow().run()
