"""Build a tool in Python, then run it in a workflow."""
from pathlib import Path

from sophios.api.python.tool_builder import (
    CommandLineTool, Input, Inputs, Output, Outputs, cwl,
)
from sophios.api.python.workflow import Step, Workflow

#: The tools that ship with Sophios, where `sophios --generate_config` put them.
ADAPTERS = Path.home() / 'wic' / 'cwl_adapters'

word_count = (
    CommandLineTool('word_count',
                    Inputs(text=Input(cwl.file, position=1)),
                    Outputs(counts=Output(cwl.file, glob='counts.txt')))
    .base_command('wc', '-w')
    .stdout('counts.txt')
    .docker('debian:bookworm-slim')
)


def workflow() -> Workflow:
    """Create a file, write a line into it, and count its words."""
    touch = Step(clt_path=ADAPTERS / 'touch.cwl')
    touch.inputs.filename = 'greeting.txt'
    append = Step(clt_path=ADAPTERS / 'append.cwl')
    append.inputs.file = touch.outputs.file
    append.inputs.str = 'Hello, Sophios'
    count = Step(word_count)
    count.inputs.text = append.outputs.file
    return Workflow([touch, append, count], 'count_words')


if __name__ == '__main__':
    word_count.write_cwl('word_count.cwl', validate=True)
    workflow().run()
