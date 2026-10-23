"""A graph: the greeting feeds two steps, and a third joins what they made."""
from pathlib import Path

from sophios.api.python.workflow import Step, Workflow


#: The tools that ship with Sophios, where `sophios --generate_config` put them.
ADAPTERS = Path.home() / 'wic' / 'cwl_adapters'


def workflow() -> Workflow:
    """Shout the greeting and count its words, side by side, then put both in one report."""
    touch = Step(clt_path=ADAPTERS / 'touch.cwl')
    touch.inputs.filename = 'greeting.txt'
    append = Step(clt_path=ADAPTERS / 'append.cwl')
    append.inputs.file = touch.outputs.file
    append.inputs.str = 'Hello, Sophios'
    shout = Step(clt_path=Path(__file__).with_name('shout.cwl'))
    shout.inputs.text = append.outputs.file
    word_count = Step(clt_path=Path(__file__).with_name('word_count.cwl'))
    word_count.inputs.text = append.outputs.file
    combine = Step(clt_path=Path(__file__).with_name('combine.cwl'))
    combine.inputs.text = shout.outputs.shouted
    combine.inputs.counts = word_count.outputs.counts
    return Workflow([touch, append, shout, word_count, combine], 'greeting_graph')


if __name__ == '__main__':
    workflow().run()
