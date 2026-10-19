"""One tool used twice: two appends, each a step under its own name."""
from pathlib import Path

from sophios.api.python.workflow import Step, Workflow


#: The tools that ship with Sophios, where `sophios --generate_config` put them.
ADAPTERS = Path.home() / 'wic' / 'cwl_adapters'


def workflow() -> Workflow:
    """Write two lines into greeting.txt with append.cwl, as the steps hello and goodbye."""
    touch = Step(clt_path=ADAPTERS / 'touch.cwl')
    touch.inputs.filename = 'greeting.txt'
    hello = Step(clt_path=ADAPTERS / 'append.cwl', step_name='hello')
    hello.inputs.str = 'Hello, Sophios'
    goodbye = Step(clt_path=ADAPTERS / 'append.cwl', step_name='goodbye')
    goodbye.inputs.str = 'Goodbye for now'
    cat = Step(clt_path=ADAPTERS / 'cat.cwl')
    return Workflow([touch, hello, goodbye, cat], 'two_lines')


if __name__ == '__main__':
    workflow().run()
