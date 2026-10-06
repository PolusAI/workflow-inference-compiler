"""Draw a compiled workflow, for `--graphviz`.

A projection made after compiling, and only when a drawing is asked for: the
compiled graph gives the steps and the edges between them, and the documents it
was compiled from their `wic: graphviz:` blocks. Nothing in the compiler reads a drawing.
"""
import graphviz

from .ir.artifacts import CompilationResult
from .ir.frontdoor import SourceBundle
from .ir.names import Names
from .ir.types import EdgeOrigin, WorkflowGraph
from .wic_types import GraphSettings


def draw(compiled: CompilationResult, bundle: SourceBundle, settings: GraphSettings,
         name: str) -> graphviz.Digraph:
    """The drawing of `compiled`, the result of compiling `bundle`.

    Args:
        compiled (CompilationResult): What compiling `bundle` returned.
        bundle (SourceBundle): The documents `compiled` was compiled from.
        settings (GraphSettings): The `--graph_*` flags.
        name (str): The drawing's name, and its title when the root document gives none.

    Returns:
        graphviz.Digraph: The drawing, ready to save as DOT.
    """
    root = graphviz.Digraph(name=name)
    # newrank='True' ranks nodes globally (rather than per-cluster), which is
    # required for rank=same constraints to work across subgraphs/clusters.
    root.attr(newrank='True')
    root.attr(bgcolor='transparent')  # Useful for making slides
    root.attr(fontcolor='black' if settings['graph_dark_theme'] else 'white')
    document = bundle.parsed.document
    drawn = dict(document.sidecar.entries).get('graphviz') if document and document.sidecar else None
    with root.subgraph(name=f'cluster_{name}') as cluster:
        cluster.attr(label=drawn.get('label', name) if isinstance(drawn, dict) else name)
        cluster.attr(color='lightblue')  # color of cluster subgraph outline
        _draw_graph(cluster, compiled.graph, Names.of(compiled.graph), settings)
    return root


def _draw_graph(drawn: graphviz.Digraph, graph: WorkflowGraph, names: Names, settings: GraphSettings) -> None:
    """Draw `graph`'s steps and edges into `drawn`, and each workflow it calls as a cluster of its own."""
    for step in graph.steps:
        label = names.step(step.id) if settings['graph_label_stepname'] else step.id.name
        drawn.node(names.qualified(step.id), label=label, shape='box', style='rounded, filled',
                   fillcolor='lightblue')
    seen: set[tuple[str, str]] = set()
    for edge in graph.edges:
        source = names.qualified(edge.source.step)
        sink = names.qualified(edge.sink.step)
        # Explicit edges are blue; an inferred one takes the theme's font
        # colour, as docs/tutorials/multistep.md says, so a reader can tell
        # what the document said from what the compiler decided.
        attrs = {'color': 'blue'}
        if edge.origin is EdgeOrigin.INFERRED:
            attrs['color'] = 'black' if settings['graph_dark_theme'] else 'white'
        if settings['graph_label_edges']:
            attrs['label'] = names.port(edge.source.port)
        elif (source, sink) in seen:
            continue
        seen.add((source, sink))
        if source != sink:
            drawn.edge(source, sink, **attrs)
    if len(graph.namespace.parts) < settings['graph_inline_depth']:
        for child in graph.children:
            cluster = graphviz.Digraph(name=f'cluster_{child.name}')
            _draw_graph(cluster, child, names, settings)
            drawn.subgraph(cluster)
