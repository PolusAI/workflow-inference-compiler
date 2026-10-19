"""An inline literal: the value is written into the workflow."""
from pathlib import Path

from sophios.api.python.workflow import Step, Workflow


#: The tools that ship with Sophios, where `sophios --generate_config` put them.
ADAPTERS = Path.home() / 'wic' / 'cwl_adapters'


def workflow() -> Workflow:
    """Touch a file whose name is given in the workflow."""
    touch = Step(clt_path=ADAPTERS / 'touch.cwl')
    touch.inputs.filename = 'empty.txt'
    return Workflow([touch], 'inline_literal')


if __name__ == '__main__':
    workflow().run()
