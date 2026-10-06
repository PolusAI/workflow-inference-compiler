"""The graph Sophios draws: what each graph setting and the CLI's `--graphviz` do to it.

These compile in memory and read the DOT source `drawing.draw` makes of the result, or
run the CLI and read the `.gv` it writes. None of them renders an image.
"""
import re
import sys
from collections.abc import Callable
from pathlib import Path
from types import SimpleNamespace

import pytest

from sophios import compiler, drawing, main, preflight
from sophios.cli import default_compilation_settings, get_args, get_dicts_for_compilation
from sophios.wic_types import Yaml

from .hermetic import bundle, subworkflow_step
from .synthetic_tools import SYNTHETIC_TOOLS

#: `mk_file` feeds `xform` of a subworkflow, which feeds `count`, which feeds a
#: subworkflow of its own: two levels of nesting under the root.
_DOCUMENT: Yaml = {'steps': [
    {'id': 'mk_file', 'in': {'name': {'wic_inline_input': 'a'}}, 'out': [{'file': {'wic_anchor': 'f'}}]},
    subworkflow_step('sub.wic', {'steps': [
        {'id': 'xform', 'in': {'file': {'wic_alias': 'f'}, 'name': {'wic_inline_input': 'b'}}},
        {'id': 'count'},
        subworkflow_step('inner.wic', {'steps': [
            {'id': 'mk_text', 'in': {'name': {'wic_inline_input': 'c'}}}]}),
    ]}),
]}


def _draw(document: Yaml = _DOCUMENT, *, label_edges: bool = False, label_stepname: bool = False,
          inline_depth: int = sys.maxsize, dark_theme: bool = False) -> str:
    """`document` compiled, then drawn with the given graph settings: the DOT source."""
    source = bundle(document, 'oracle', SYNTHETIC_TOOLS)
    compiled = compiler.compile_source(source, default_compilation_settings(), relative_run_path=True, testing=True)
    _options, settings = get_dicts_for_compilation(get_args())
    settings['graph_label_edges'] = label_edges
    settings['graph_label_stepname'] = label_stepname
    settings['graph_inline_depth'] = inline_depth
    settings['graph_dark_theme'] = dark_theme
    dot: str = drawing.draw(compiled, source, settings, 'oracle').source
    return dot


#: A node or edge name in DOT: quoted, or a bare word.
_ID = r'"[^"]*"|\w+'
#: A node or edge statement: its node or its two ends, and its attributes.
_STATEMENT = re.compile(rf'^\s*({_ID})(?: -> ({_ID}))? \[(.*)\]$', re.MULTILINE)
_ATTRIBUTE = re.compile(r'(\w+)=("(?:[^"\\]|\\.)*"|[^\s"]+)')


def _unquoted(text: str) -> str:
    return text[1:-1] if text.startswith('"') else text


def _statements(dot: str) -> list[tuple[str, str, dict[str, str]]]:
    """Each node (with an empty sink) and each edge `dot` draws, with its attributes."""
    return [(_unquoted(name), _unquoted(sink), {key: _unquoted(value) for key, value in _ATTRIBUTE.findall(attrs)})
            for name, sink, attrs in _STATEMENT.findall(dot)]


def _nodes(dot: str) -> dict[str, dict[str, str]]:
    return {name: attrs for name, sink, attrs in _statements(dot) if not sink}


def _edges(dot: str) -> list[tuple[str, str, dict[str, str]]]:
    return [statement for statement in _statements(dot) if statement[1]]


@pytest.mark.fast
@pytest.mark.parametrize(('labelled', 'expected'), [(True, 'file'), (False, None)])
def test_graph_label_edges_names_the_output_each_edge_carries(labelled: bool, expected: str | None) -> None:
    """With `--graph_label_edges` an edge is labelled with the output it carries."""
    dot = _draw(label_edges=labelled)
    edges = _edges(dot)
    assert len(edges) == 2
    assert [attrs.get('label') for _source, _sink, attrs in edges] == [expected] * 2
    assert ('label=file' in dot) is labelled


@pytest.mark.fast
@pytest.mark.parametrize(('qualified', 'labels'), [
    (False, {'mk_file', 'sub.wic', 'xform', 'count', 'inner.wic', 'mk_text'}),
    (True, {'oracle__step__1__mk_file', 'oracle__step__2__sub.wic', 'sub__step__1__xform',
            'sub__step__2__count', 'sub__step__3__inner.wic', 'inner__step__1__mk_text'}),
])
def test_graph_label_stepname_labels_a_node_with_the_generated_step_name(qualified: bool, labels: set[str]) -> None:
    """A node is labelled with the id the document wrote; with `--graph_label_stepname` it is the generated
    step name instead, which a subworkflow's steps carry relative to the subworkflow."""
    assert {attrs['label'] for attrs in _nodes(_draw(label_stepname=qualified)).values()} == labels


@pytest.mark.fast
@pytest.mark.parametrize(('depth', 'clusters'), [
    (0, {'cluster_oracle'}),
    (1, {'cluster_oracle', 'cluster_sub'}),
    (2, {'cluster_oracle', 'cluster_sub', 'cluster_inner'}),
])
def test_graph_inline_depth_draws_that_many_levels_of_subworkflow(depth: int, clusters: set[str]) -> None:
    """`--graph_inline_depth N` draws the subworkflows N levels down as clusters and leaves the steps of the
    deeper ones out of the drawing."""
    dot = _draw(inline_depth=depth)
    assert set(re.findall(r'subgraph (cluster_\w+)', dot)) == clusters
    drawn = _nodes(dot)
    assert ('oracle__step__2__sub.wic___sub__step__2__count' in drawn) is (depth >= 1)
    assert ('oracle__step__2__sub.wic___sub__step__3__inner.wic___inner__step__1__mk_text' in drawn) is (depth >= 2)


@pytest.mark.fast
@pytest.mark.parametrize(('dark_theme', 'font_colour'), [(False, 'white'), (True, 'black')])
def test_the_graph_draws_inferred_edges_in_the_font_colour_and_the_rest_in_blue(
        dark_theme: bool, font_colour: str) -> None:
    """An edge the document wrote is blue, whether it joins two steps of one
    document or reaches into a subworkflow; an edge the compiler inferred takes
    the theme's font colour, which is black on a dark theme and white otherwise.
    """
    document: Yaml = {'steps': [
        {'id': 'mk_file', 'in': {'name': {'wic_inline_input': 'a'}}, 'out': [{'file': {'wic_anchor': 'f'}}]},
        {'id': 'xform', 'in': {'file': {'wic_alias': 'f'}, 'name': {'wic_inline_input': 'b'}}},
        subworkflow_step('sub.wic', {'steps': [
            {'id': 'xform', 'in': {'file': {'wic_alias': 'f'}, 'name': {'wic_inline_input': 'c'}}},
            {'id': 'count'},                                       # file inferred from xform
        ]}),
    ]}
    colours = {(source, sink): attrs.get('color')
               for source, sink, attrs in _edges(_draw(document, dark_theme=dark_theme))}
    assert colours == {
        ('oracle__step__1__mk_file', 'oracle__step__2__xform'): 'blue',
        ('oracle__step__1__mk_file', 'oracle__step__3__sub.wic___sub__step__1__xform'): 'blue',
        ('oracle__step__3__sub.wic___sub__step__1__xform',
         'oracle__step__3__sub.wic___sub__step__2__count'): font_colour,
    }


@pytest.fixture(name='cli_graph')
def _cli_graph(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> Callable[..., str]:
    """Run the CLI on a workflow in `tmp_path` and return the `.gv` it wrote.

    `dot` is reported present and the commands that would run it are swallowed: what is under
    test is the DOT source Sophios writes, not the image Graphviz makes of it.
    """
    monkeypatch.setattr(preflight, 'prepare', lambda *_a, **_k: None)
    monkeypatch.setattr(main.shutil, 'which', lambda _name: '/usr/bin/dot')
    monkeypatch.setattr(main, 'sub', SimpleNamespace(run=lambda *_a, **_k: None))
    monkeypatch.chdir(tmp_path)

    def run(source: str, *flags: str) -> str:
        workflow = tmp_path / 'w.wic'
        workflow.write_text(source, encoding='utf-8')
        monkeypatch.setattr('sys.argv', ['sophios', '--yaml', str(workflow), '--generate_cwl_workflow',
                                         '--graphviz', *flags])
        main.main()
        return (tmp_path / 'autogenerated' / 'w.wic.gv').read_text(encoding='utf-8')
    return run


_TWO_STEPS = """\
steps:
- id: touch
  in:
    filename: !ii empty.txt
- id: cat
"""


@pytest.mark.fast
def test_graphviz_writes_the_graph_of_the_compiled_workflow(cli_graph: Callable[..., str]) -> None:
    """`--graphviz` writes `autogenerated/<stem>.wic.gv` with a node per step and an edge per dependency."""
    drawn = cli_graph(_TWO_STEPS)
    assert re.findall(r'^\s*(\w+) \[label=(\w+) ', drawn, re.MULTILINE) == [
        ('w__step__1__touch', 'touch'), ('w__step__2__cat', 'cat')]
    assert re.findall(r'^\s*(\w+) -> (\w+) \[', drawn, re.MULTILINE) == [('w__step__1__touch', 'w__step__2__cat')]


@pytest.mark.fast
@pytest.mark.parametrize(('flags', 'font'), [((), 'white'), (('--graph_dark_theme',), 'black')])
def test_graph_dark_theme_sets_the_font_of_the_drawing(
        cli_graph: Callable[..., str], flags: tuple[str, ...], font: str) -> None:
    """The drawing's font and its inferred edges take one colour, white, or black with `--graph_dark_theme`."""
    drawn = cli_graph(_TWO_STEPS, *flags)
    assert f'fontcolor={font}\n' in drawn
    assert f'-> w__step__2__cat [color={font}]' in drawn


@pytest.mark.fast
def test_the_root_graphviz_label_titles_the_drawing(cli_graph: Callable[..., str]) -> None:
    """`wic: graphviz: label:` on the workflow is the title of the cluster that draws it."""
    drawn = cli_graph('wic:\n  graphviz:\n    label: My Pipeline\n' + _TWO_STEPS)
    assert 'label="My Pipeline"' in drawn
    assert 'label="My Pipeline"' not in cli_graph(_TWO_STEPS)
