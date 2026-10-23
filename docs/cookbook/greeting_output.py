"""The greeting file, returned under a name of its own: a workflow output."""
from pathlib import Path

from sophios.api.python.workflow import Step, Workflow


#: The tools that ship with Sophios, where `sophios --generate_config` put them.
ADAPTERS = Path.home() / 'wic' / 'cwl_adapters'


def workflow() -> Workflow:
    """The greeting workflow, with append's file named ``greeting``."""
    touch = Step(clt_path=ADAPTERS / 'touch.cwl')
    touch.inputs.filename = 'greeting.txt'
    append = Step(clt_path=ADAPTERS / 'append.cwl')
    append.inputs.file = touch.outputs.file
    append.inputs.str = 'Hello, Sophios'
    cat = Step(clt_path=ADAPTERS / 'cat.cwl')
    cat.inputs.file = append.outputs.file
    flow = Workflow([touch, append, cat], 'greeting_output')
    flow.outputs.greeting = append.outputs.file
    return flow


if __name__ == '__main__':
    workflow().run()
