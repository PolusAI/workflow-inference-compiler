"""The greeting with its connections left out: edge inference fills them in."""
from pathlib import Path

from sophios.api.python.workflow import Step, Workflow


#: The tools that ship with Sophios, where `sophios --generate_config` put them.
ADAPTERS = Path.home() / 'wic' / 'cwl_adapters'


def workflow() -> Workflow:
    """The same three steps as greeting.py, with no input connected by hand."""
    touch = Step(clt_path=ADAPTERS / 'touch.cwl')
    touch.inputs.filename = 'greeting.txt'
    append = Step(clt_path=ADAPTERS / 'append.cwl')
    append.inputs.str = 'Hello, Sophios'
    cat = Step(clt_path=ADAPTERS / 'cat.cwl')
    return Workflow([touch, append, cat], 'greeting_inferred')


if __name__ == '__main__':
    workflow().run()
