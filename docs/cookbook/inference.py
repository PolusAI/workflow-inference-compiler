"""Inference: an input left unbound is connected by the compiler."""
from pathlib import Path

from sophios.api.python.workflow import Step, Workflow


#: The tools that ship with Sophios, where `sophios --generate_config` put them.
ADAPTERS = Path.home() / 'wic' / 'cwl_adapters'


def workflow() -> Workflow:
    """Print a file; which file is left to inference."""
    touch = Step(clt_path=ADAPTERS / 'touch.cwl')
    touch.inputs.filename = 'empty.txt'
    echo = Step(clt_path=ADAPTERS / 'echo.cwl')
    echo.inputs.message = 'Hello'
    cat = Step(clt_path=ADAPTERS / 'cat.cwl')
    return Workflow([touch, echo, cat], 'inference')


if __name__ == '__main__':
    workflow().run()
