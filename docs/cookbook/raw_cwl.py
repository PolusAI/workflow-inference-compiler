"""Raw CWL: the name is written into the CWL as it is."""
from pathlib import Path

from sophios.api.python.workflow import Step, Workflow


#: The tools that ship with Sophios, where `sophios --generate_config` put them.
ADAPTERS = Path.home() / 'wic' / 'cwl_adapters'


def workflow() -> Workflow:
    """Echo the workflow input ``greeting``."""
    echo = Step(clt_path=ADAPTERS / 'echo.cwl')
    flow = Workflow([echo], 'raw_cwl')
    echo.inputs.message = flow.inputs.greeting
    return flow


if __name__ == '__main__':
    workflow().run()
