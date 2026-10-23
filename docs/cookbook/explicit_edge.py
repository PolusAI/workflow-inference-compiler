"""An explicit edge: one step's output feeds a later step's input."""
from pathlib import Path

from sophios.api.python.workflow import Step, Workflow


#: The tools that ship with Sophios, where `sophios --generate_config` put them.
ADAPTERS = Path.home() / 'wic' / 'cwl_adapters'


def workflow() -> Workflow:
    """Touch a file, then append a line to that file."""
    touch = Step(clt_path=ADAPTERS / 'touch.cwl')
    touch.inputs.filename = 'hello.txt'
    append = Step(clt_path=ADAPTERS / 'append.cwl')
    append.inputs.file = touch.outputs.file
    append.inputs.str = 'Hello'
    return Workflow([touch, append], 'explicit_edge')


if __name__ == '__main__':
    workflow().run()
