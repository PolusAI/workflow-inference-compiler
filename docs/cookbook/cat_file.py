"""A File literal: the workflow names a file that exists before it runs."""
from pathlib import Path

from sophios.api.python.workflow import Step, Workflow


#: The tools that ship with Sophios, where `sophios --generate_config` put them.
ADAPTERS = Path.home() / 'wic' / 'cwl_adapters'


def workflow() -> Workflow:
    """Print data.txt, read from the working directory."""
    cat = Step(clt_path=ADAPTERS / 'cat.cwl')
    cat.inputs.file = {'class': 'File', 'location': 'data.txt'}
    return Workflow([cat], 'cat_file')


if __name__ == '__main__':
    workflow().run()
