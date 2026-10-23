"""A renaming: the abstract step's input_abstract should read output_summary."""
from pathlib import Path

from sophios.api.python.workflow import Step, Workflow


#: The tools that ship with Sophios, where `sophios --generate_config` put them.
ADAPTERS = Path.home() / 'wic' / 'cwl_adapters'


def workflow() -> Workflow:
    """Write a text, summarize it, and write an abstract of it, every connection left to inference."""
    echo = Step(clt_path=ADAPTERS / 'echo.cwl')
    echo.inputs.message = 'The quick brown fox jumps over the lazy dog'
    summarize = Step(clt_path=Path(__file__).with_name('summarize.cwl'))
    abstract = Step(clt_path=Path(__file__).with_name('abstract.cwl'))
    return Workflow([echo, summarize, abstract], 'naming_renamed')


if __name__ == '__main__':
    workflow().run()
