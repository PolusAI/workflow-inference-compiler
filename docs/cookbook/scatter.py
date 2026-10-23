"""Scatter: run a step once per combination of array elements."""
from pathlib import Path

from sophios.api.python.workflow import Step, Workflow


#: The tools that ship with Sophios, where `sophios --generate_config` put them.
ADAPTERS = Path.home() / 'wic' / 'cwl_adapters'


def workflow() -> Workflow:
    """Echo every pairing of two lists of words."""
    echo_3 = Step(clt_path=ADAPTERS / 'echo_3.cwl')
    echo_3.inputs.message1 = ['Hello', 'Goodbye']
    echo_3.inputs.message2 = ['sun', 'moon']
    echo_3.inputs.message3 = 'and stars'
    echo_3.scatter_on(echo_3.inputs.message1, echo_3.inputs.message2, method='flat_crossproduct')
    return Workflow([echo_3], 'scatter')


if __name__ == '__main__':
    workflow().run()
