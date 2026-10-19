"""Scatter: the greeting pipeline, run once for each of three people."""
from pathlib import Path

from sophios.api.python.workflow import Step, Workflow


#: The tools that ship with Sophios, where `sophios --generate_config` put them.
ADAPTERS = Path.home() / 'wic' / 'cwl_adapters'


def workflow() -> Workflow:
    """Make, write and read one greeting file per name, the lists walked side by side."""
    touch = Step(clt_path=ADAPTERS / 'touch.cwl')
    touch.inputs.filename = ['alice.txt', 'bob.txt', 'carol.txt']
    touch.scatter_on(touch.inputs.filename)
    append = Step(clt_path=ADAPTERS / 'append.cwl')
    append.inputs.file = touch.outputs.file
    append.inputs.str = ['Hello, Alice', 'Hello, Bob', 'Hello, Carol']
    append.scatter_on(append.inputs.file, append.inputs.str)
    cat = Step(clt_path=ADAPTERS / 'cat.cwl')
    cat.inputs.file = append.outputs.file
    cat.scatter_on(cat.inputs.file)
    return Workflow([touch, append, cat], 'greet_many')


if __name__ == '__main__':
    workflow().run()
