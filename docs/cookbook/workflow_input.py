"""A workflow input: the value is given when the workflow runs."""
from pathlib import Path

from sophios.api.python.workflow import Step, Workflow


#: The tools that ship with Sophios, where `sophios --generate_config` put them.
ADAPTERS = Path.home() / 'wic' / 'cwl_adapters'


def workflow() -> Workflow:
    """Touch a file whose name is the workflow input ``filename``."""
    touch = Step(clt_path=ADAPTERS / 'touch.cwl')
    flow = Workflow([touch], 'workflow_input')
    touch.inputs.filename = flow.inputs.filename
    return flow


if __name__ == '__main__':
    workflow().run()
