"""`Workflow.draw` returns the drawing `--graphviz` makes, on a background of its own, from the same compile."""
import re
from pathlib import Path

import graphviz
import pytest

from sophios.api.python.workflow import Inference, Step, Workflow

ADAPTERS = Path(__file__).resolve().parents[2] / "cwl_adapters"


def _workflow() -> Workflow:
    """`touch` feeds `append` by an edge written here, and `append` feeds `cat` by an edge inference adds."""
    touch = Step(clt_path=ADAPTERS / "touch.cwl")
    touch.inputs.filename = "empty.txt"
    append = Step(clt_path=ADAPTERS / "append.cwl")
    append.inputs.str = "Hello"
    append.inputs.file = touch.outputs.file
    cat = Step(clt_path=ADAPTERS / "cat.cwl")
    return Workflow([touch, append, cat], "wf")


def _box_labels(drawing: graphviz.Digraph) -> list[str]:
    """The label of each box, in the order the drawing declares the boxes."""
    return re.findall(r'\[label="?([\w@~]+)"? ', drawing.source)


def _edge_colours(drawing: graphviz.Digraph) -> dict[tuple[str, str], str]:
    """The colour of each edge, by the labels of the steps it joins."""
    edges = re.findall(r'^\s*"?(\w+)"? -> "?(\w+)"? \[color=(\w+)\]$',
                       drawing.source, re.MULTILINE)
    return {(source, sink): colour for source, sink, colour in edges}


@pytest.mark.fast
def test_draw_paints_a_white_page_and_colours_the_edges_by_who_made_them() -> None:
    """By default the drawing is on an opaque white background, in black, and an edge written
    in the workflow is blue where one inference added is black."""
    drawing = _workflow().draw()
    assert isinstance(drawing, graphviz.Digraph)
    assert 'bgcolor=white' in drawing.source
    assert 'fontcolor=black' in drawing.source
    assert _edge_colours(drawing) == {('touch', 'append'): 'blue', ('append', 'cat'): 'black'}
    assert 'transparent' not in drawing.source


@pytest.mark.fast
def test_draw_dark_paints_a_dark_page_in_white_with_a_lighter_blue() -> None:
    """`dark=True` is an opaque near-black background, white text and outlines, a white inferred edge
    and a lighter blue for the edge written in the workflow."""
    drawing = _workflow().draw(dark=True)
    assert 'bgcolor="#1e1e1e"' in drawing.source
    assert 'fontcolor=white' in drawing.source
    assert re.findall(r'\[label=\w+ color=(\w+) fontcolor=(\w+) shape=box', drawing.source) == [
        ('white', 'white')] * 3
    assert _edge_colours(drawing) == {('touch', 'append'): 'lightskyblue', ('append', 'cat'): 'white'}


@pytest.mark.fast
def test_draw_infers_as_compile_does() -> None:
    """With `Inference(enabled=False)` the edge inference adds is not drawn; the one written is."""
    assert _edge_colours(_workflow().draw()) == {('touch', 'append'): 'blue', ('append', 'cat'): 'black'}
    assert _edge_colours(_workflow().draw(inference=Inference(enabled=False))) == {('touch', 'append'): 'blue'}


@pytest.mark.fast
def test_draw_labels_edges_and_boxes_as_the_flags_do() -> None:
    """`label_edges` names each edge's port, and `label_stepname` names each box by its name in the
    compiled CWL, which differs from the step's id when a workflow output has that id: `touch` is
    compiled as `touch@1`."""
    assert 'label=file' in _workflow().draw(label_edges=True).source
    flow = _workflow()
    flow.outputs.touch = flow.steps[0].outputs.file
    assert _box_labels(flow.draw()) == ['touch', 'append', 'cat']
    assert _box_labels(flow.draw(label_stepname=True)) == ['touch@1', 'append', 'cat']
