import argparse
import copy
from pathlib import Path
from typing import List, Tuple

import graphviz
import yaml
from typing import Any, Dict

from . import auto_gen_header
from .wic_types import GraphData, Json, Namespaces, NodeData, Yaml, RoseTree, YamlForest, GraphReps, Tool, Tools, YamlTree

# Use white for dark backgrounds, black for light backgrounds
font_edge_color = 'white'

def read_lines_pairs(filename: Path) -> List[Tuple[str, str]]:
    lines = []
    for line in open(filename, 'r').readlines():
        if line.strip() == '':  # Skip blank lines
            continue
        if line.startswith('#'):  # Skip comment lines
            continue
        ls = line.split()
        if not len(ls) == 2:
            print(line)
            raise Exception("Error! Line must contain exactly two entries!")
        lines.append((ls[0], ls[1]))
    return lines


def step_name_str(yaml_stem: str, i: int, step_key: str) -> str:
    """Returns a string which uniquely and hierarchically identifies a step in a workflow

    Args:
        yaml_stem (str): The name of the workflow (filepath stem)
        i (int): The (zero-based) step number
        step_key (str): The name of the step (used as a dict key)

    Returns:
        str: The parameters (and the word 'step') joined together with double underscores
    """
    # Use double underscore so we can '__'.split() below.
    # (This should work as long as yaml_stem and step_key do not contain __)
    return f'{yaml_stem}__step__{i+1}__{step_key}'


def parse_step_name_str(step_name_str: str) -> Tuple[str, int, str]:
    """The inverse function to step_name_str()

    Args:
        step_name_str (str): A string of the same form as returned by step_name_str()

    Raises:
        Exception: If the argument is not of the same form as returned by step_name_str()

    Returns:
        Tuple[str, int, str]: The parameters used to create step_name_str
    """
    vals = step_name_str.split('__') # double underscore
    if not len(vals) == 4:
        raise Exception(f"Error! {step_name_str} is not of the format \n"
                        + '{yaml_stem}__step__{i+1}__{step_key}\n'
                        + 'yaml_stem and step_key should not contain any double underscores.')
    try:
        i = int(vals[2])
    except:
        raise Exception(f"Error! {step_name_str} is not of the format \n"
                        + '{yaml_stem}__step__{i+1}__{step_key}')
    return (vals[0], i-1, vals[3])


def add_yamldict_keyval(steps_i: Yaml, step_key: str, in_out: str, keyval: Yaml) -> Yaml:
    """Convenience function used to (mutably) merge two Yaml dicts.

    Args:
        steps_i (Yaml): A partially-completed Yaml dict representing a step in a CWL workflow
        step_key (str): The name of the step in a CWL workflow
        in_out (str): Either the string 'in' or the string 'out'
        keyval (Yaml): A Yaml dict with additional details to be merged into the first Yaml dict

    Returns:
        Yaml: The first Yaml dict with the second Yaml dict merged into it.
    """
    # TODO: Check whether we can just use deepmerge.merge()
    if steps_i[step_key]:
        if in_out in steps_i[step_key]:
            new_keys = dict(list(steps_i[step_key][in_out].items()) + list(keyval.items()))
            new_keyvals = dict([(k, v) if k != in_out else (k, new_keys) for k, v in steps_i[step_key].items()])
        else:
            new_keys = keyval
            new_keyvals = dict(list(steps_i[step_key].items()) + [(in_out, new_keys)])
        steps_i[step_key].update(new_keyvals)
    else:
        steps_i = {step_key: {in_out: keyval}}
    return steps_i

def add_yamldict_keyval_in(steps_i: Yaml, step_key: str, keyval: Yaml) -> Yaml:
    return add_yamldict_keyval(steps_i, step_key, 'in', keyval)

def add_yamldict_keyval_out(steps_i: Yaml, step_key: str, keyval: List[str]) -> Yaml:
    return add_yamldict_keyval(steps_i, step_key, 'out', keyval) # type: ignore


def add_graph_edge(args: argparse.Namespace, graph: GraphReps, nss1: Namespaces, nss2: Namespaces, label: str, color: str = font_edge_color) -> None:
    """Adds edges to (all of) our graph representations, with the ability to
    collapse all nodes below a given depth to a single node.

    This function utilizes the fact that nodes have been carefully designed to
    have unique, hierarchical names. If we want to hide all of the details
    below a given depth, we can simply truncate each of the namespaces!
    (and do the same when creating the nodes)

    Args:
        args (argparse.Namespace): The command line arguments
        graph (GraphReps): A tuple of a GraphViz DiGraph and a networkx DiGraph
        nss1 (Namespaces): The namespaces associated with the first node
        nss2 (Namespaces): The namespaces associated with the second node
        label (str): The edge label
        color (str, optional): The edge color Defaults to font_edge_color
    """
    nss1 = nss1[:(1 + args.graph_inline_depth)]
    edge_node1 = '___'.join(nss1)
    nss2 = nss2[:(1 + args.graph_inline_depth)]
    edge_node2 = '___'.join(nss2)
    graph_gv = graph.graphviz
    graph_nx = graph.networkx
    graphdata = graph.graphdata
    # Hide internal self-edges
    if not edge_node1 == edge_node2:
        attrs = {'color': color}
        if args.graph_label_edges:
            attrs['label'] = label

        graph_gv.edge(edge_node1, edge_node2, **attrs)
    graph_nx.add_edge(edge_node1, edge_node2)
    graphdata.edges.append((edge_node1, edge_node2, attrs))


def flatten_graphdata(graphdata: GraphData, parent: str = '') -> GraphData:
    subgraphs = [flatten_graphdata(subgraph, str(graphdata.name)) for subgraph in graphdata.subgraphs]

    # NOTE: Even though all of the following default list arguments are [],
    # you MUST explicitly supply the empty lists!!! Otherwise, after
    # instantiation, the lists will contain values from previous instances!!!
    # This shallow copy causes an infinite loop because as we copy nodes,
    # they end up getting appended to the original lists!
    # This makes absolutely no sense. Since the lists are defined at the
    # instance level (NOT the class level), there should be zero sharing!
    gd = GraphData(str(graphdata.name), [], [], [], []) # This is fine
    #gd = GraphData(str(graphdata.name)) # This is NOT fine!
    # i.e. The following statement will NOT print zeros (in the second case)!
    #print(gd.name, len(gd.nodes), len(gd.edges), len(gd.subgraphs), len(gd.ranksame))

    for subgraph in subgraphs:
        # We need to add a placeholder node for each subgraph first
        attrs = {} if parent == '' else {'parent': parent}
        # NOTE: This does not yet work with args.graph_inline_depth
        gd.nodes.append((subgraph.name, attrs))

    for subgraph in subgraphs:
        # Then we can add the nodes and edges from the subgraphs.
        # (Otherwise, cytoscape won't render the subgraphs correctly.)
        for (subnode, subattrs) in subgraph.nodes:
            gd.nodes.append((subnode, subattrs))
        for (subnode1, subnode2, subattrs) in subgraph.edges:
            gd.edges.append((subnode1, subnode2, subattrs))

    # Finally, add the nodes and edges from the current graph
    for (node, attrs) in graphdata.nodes:
        attrs['parent'] = graphdata.name
        gd.nodes.append((node, attrs))
    for (node1, node2, attrs) in graphdata.edges:
        gd.edges.append((node1, node2, attrs))

    return gd


def graphdata_to_cytoscape(graphdata: GraphData) -> Json:
    nodes = []
    for (node, attrs) in list(graphdata.nodes):
        nodes.append({'data': {'id': node, **attrs}})
    edges = []
    for (node1, node2, attrs) in list(graphdata.edges):
        edges.append({'data': {'source': node1, 'target': node2, **attrs}})
    return {'nodes': nodes, 'edges': edges}


def partition_by_lowest_common_ancestor(nss1: Namespaces, nss2: Namespaces) -> Tuple[Namespaces, Namespaces]:
    """See https://en.wikipedia.org/wiki/Lowest_common_ancestor

    Args:
        nss1 (Namespaces): The namespaces associated with the first node
        nss2 (Namespaces): The namespaces associated with the second node

    Returns:
        Tuple[Namespaces, Namespaces]: nss1, partitioned by lowest common ancestor
    """
    # Only partition nss1; if you want to partition nss1
    # just switch the arguments at the call site.
    if nss1 == [] or nss2 == []:
        return ([], nss1) # Base case
    if nss1[0] == nss2[0]: # Keep going
        (nss1_heads, nss1_tails) = partition_by_lowest_common_ancestor(nss1[1:], nss2[1:])
        return ([nss1[0]] + nss1_heads, nss1_tails)
    else:
        return ([], nss1)


def get_steps_keys(steps: List[Yaml]) -> List[str]:
    """Returns the name (dict key) of each step in the given CWL workflow

    Args:
        steps (List[Yaml]): The steps: tag of a CWL workflow

    Returns:
        List[str]: The name of each step in the given CWL workflow
    """
    # Get the dictionary key (i.e. the name) of each step.
    steps_keys = []
    for step in steps:
        steps_keys += list(step)
    #print(steps_keys)
    return steps_keys


def extract_backend(yaml_tree: Yaml, wic: Yaml, yaml_path: Path) -> Tuple[str, Yaml]:
    """Chooses a specific backend for a given CWL workflow step.

    The backends should be thought of as either 'exactly' identical, or at
    least the same high-level protocol but implemented with a different algorithm.

    Args:
        yaml_tree (Yaml): A Yaml AST dict with sub-dicts for each backend.
        yaml_path (Path): The filepath of yaml_tree, only used for error reporting.

    Raises:
        Exception: If the steps: and/or backend: tags are not present.

    Returns:
        Tuple[str, Yaml]: The Yaml AST dict of the chosen backend.
    """
    yaml_tree_copy = copy.deepcopy(yaml_tree)
    backend = ''
    if 'backends' in wic:
        if 'default_backend' in wic:
            backend = wic['default_backend']
        if 'backend' in wic:
            backend = wic['backend']
        if backend == '':
            raise Exception(f'Error! No backend in {yaml_path}!')
        if backend not in wic['backends']:
            raise Exception(f'Error! No steps for backend {backend} in {yaml_path}!')
        steps = wic['backends'][backend]['steps']
        yaml_tree_copy.update({'steps': steps})
    elif 'steps' in yaml_tree_copy:
        pass # steps = yaml_tree_copy['steps']
    else:
        raise Exception(f'Error! No backends and/or steps in {yaml_path}!')
    return (backend, yaml_tree_copy)


def inline_sub_steps(yaml_path: Path, tools: Tools, yml_paths: Dict[str, Path]) -> List[Yaml]:
    """Recursively inlines the contents of ALL of the yml sub-workflows. (deprecated)
    
    This function is deprecated and will soon be replaced with a better implementation.

    Args:
        yaml_path (Path): The filepath of the yml workflow.
        tools (Tools): The CWL CommandLineTool definitions found using get_tools_cwl()
        yml_paths (Dict[str, Path]): The yml workflow definitions found using get_yml_paths()

    Returns:
        List[Yaml]: The recursively inlined contents of the given yml workflow.
    """
    # Load the high-level yaml workflow file.
    with open(Path(yaml_path), 'r') as y:
        yaml_tree: Yaml = yaml.safe_load(y.read())

    wic = yaml_tree.get('wic', {})
    (back_name_, yaml_tree) = extract_backend(yaml_tree, wic, yaml_path)
    steps = yaml_tree['steps']

    # Get the dictionary key (i.e. the name) of each step.
    steps_keys = []
    for step in steps:
        steps_keys += list(step)

    subkeys = [key for key in steps_keys if key not in tools]

    steps_all = []
    for i, step_key in enumerate(steps_keys):
        if step_key in subkeys:
            path = yml_paths[Path(step_key).stem]
            steps_i = inline_sub_steps(path, tools, yml_paths)
        else:
            steps_i = [steps[i]]
        steps_all.append(steps_i)

    steps_all_flattened = [step for steps in steps_all for step in steps]
    return steps_all_flattened


def flatten(lists: List[List[Any]]) -> List[Any]:
    """Concatenates a list of lists into a single list.

    Args:
        lists (List[List[Any]]): A list of lists

    Returns:
        List[Any]: A single list
    """
    return [x for lst in lists for x in lst]


def flatten_rose_tree(rose_tree: RoseTree) -> List[Any]:
    """Flattens the data contained in the Rose Tree into a List

    Args:
        rose_tree (RoseTree): A Rose Tree

    Returns:
        List[Any]: The list of data associated with each node in the RoseTree
    """
    sub_rose_trees = [flatten_rose_tree(r) for r in rose_tree.sub_trees]
    return [rose_tree.data] + flatten(sub_rose_trees)


def print_forest(forest: YamlForest) -> None:
    print(forest.yaml_tree.name)
    print(yaml.dump(forest.yaml_tree.yml))
    print(yaml.dump(forest.sub_forests))


def flatten_forest(forest: YamlForest) -> List[YamlForest]:
    """Flattens the sub-trees encountered while traversing an AST

    Args:
        forest (YamlForest): The yaml AST forest to be flattened

    Raises:
        Exception: If backend: tags are missing.

    Returns:
        List[YamlForest]: The flattened forest
    """
    #print_forest(forest)
    if forest == {}:
        return []
    yaml_tree = forest.yaml_tree.yml
    wic = {'wic': yaml_tree.get('wic', {})}

    if 'backends' in wic['wic']:
        #print_forest(forest)
        back_name = ''
        if 'default_backend' in wic['wic']:
            back_name = wic['wic']['default_backend']
        if 'backend' in wic['wic']:
            back_name = wic['wic']['backend']
        if back_name == '':
            print_forest(forest)
            raise Exception(f'Error! No backend in yaml forest!\n')
        yaml_tree_back: YamlTree = forest.sub_forests[back_name].yaml_tree
        step_1 = yaml_tree_back.yml['steps'][0]
        step_name_1 = list(step_1.keys())[0]
        if Path(step_name_1).suffix == '.yml':
            # Choose a specific backend
            return flatten_forest(forest.sub_forests[back_name])
        else:
            return [forest]

    forests = list(forest.sub_forests.values())
    sub_forests = [flatten_forest(f) for f in forests]
    # Use depth first search flattening to match flatten_rose_tree()
    #bfs = forests + flatten(sub_forests)
    dfs_lists = [[f] + fs for f, fs in zip(forests, sub_forests)]
    dfs = flatten(dfs_lists)
    return dfs


def write_to_disk(rose_tree: RoseTree, path: Path, relative_run_path: bool) -> None:
    """Writes the compiled CWL files and their associated yml inputs files to disk.

    NOTE: Only the yml input file associated with the root workflow is
    guaranteed to have all inputs. In other words, subworkflows will all have
    valid CWL files, but may not be executable due to 'missing' inputs.

    Args:
        rose_tree (RoseTree): The data associated with compiled subworkflows
        path (Path): The directory in which to write the files
        relative_run_path (bool): Controls whether to use subdirectories or just one directory.
    """
    node_data: NodeData = rose_tree.data
    namespaces = node_data.namespaces
    yaml_stem = node_data.name
    cwl_tree = node_data.compiled_cwl
    yaml_inputs = node_data.workflow_inputs_file

    path.mkdir(parents=True, exist_ok=True)
    if relative_run_path:
        filename_cwl = f'{yaml_stem}.cwl'
        filename_yml = f'{yaml_stem}_inputs.yml'
    else:
        filename_cwl = '___'.join(namespaces + [f'{yaml_stem}.cwl'])
        filename_yml = '___'.join(namespaces + [f'{yaml_stem}_inputs.yml'])

    # Dump the compiled CWL file contents to disk.
    # Use sort_keys=False to preserve the order of the steps.
    yaml_content = yaml.dump(cwl_tree, sort_keys=False, line_break='\n', indent=2)
    with open(path / filename_cwl, 'w') as w:
        w.write('#!/usr/bin/env cwl-runner\n')
        w.write(auto_gen_header)
        w.write(''.join(yaml_content))

    yaml_content = yaml.dump(yaml_inputs, sort_keys=False, line_break='\n', indent=2)
    with open(path / filename_yml, 'w') as inp:
        inp.write(auto_gen_header)
        inp.write(yaml_content)

    for sub_rose_tree in rose_tree.sub_trees:
        subpath = path
        if relative_run_path:
            sub_node_data: NodeData = sub_rose_tree.data
            sub_step_name = sub_node_data.namespaces[-1]
            subpath = path / sub_step_name
        write_to_disk(sub_rose_tree, subpath, relative_run_path)


def recursively_delete_dict_key(key: str, obj: Any) -> Any:
    """Recursively deletes any dict entries with the given key.

    Args:
        key (str): The key to be deleted
        obj (Any): The object from which to delete key.

    Returns:
        Any: _description_
    """
    if isinstance(obj, List):
        return [recursively_delete_dict_key(key, x) for x in obj]
    elif isinstance(obj, Dict):
        new_dict = {}
        for key_ in obj.keys():
            if not (key_ == key): # i.e. effectively delete key
                new_dict[key_] = recursively_delete_dict_key(key, obj[key_])
        return new_dict
    else:
        return obj


def make_tool_DAG(tool_stem: str, tool: Tool) -> None:
    """Uses the `dot` executable from the graphviz package to make a Directed
    Acyclic Graph corresponding to the given CWL CommandLineTool

    Args:
        tool_stem (str): The name of the Tool
        tool (Tool): The CWL ComandLineTool
    """
    (tool_path, tool_cwl) = tool
    yaml_path = f'autogenerated/DAG/{tool_path}'
    Path(yaml_path).parent.mkdir(parents=True, exist_ok=True)
    graph = graphviz.Digraph(name=yaml_path)
    graph.attr(bgcolor="transparent") # Useful for making slides
    graph.attr(fontcolor=font_edge_color)
    graph.attr(rankdir='LR')
    graph.node(tool_stem, shape='box', style='rounded, filled', fillcolor='lightblue')
    for input in tool_cwl['inputs']:
        input_initial_ns = input.split('___')[0].split('__')[0]
        input_no_initial_ns = input.replace(f'{input_initial_ns}__', '')
        # Hide optional inputs that could be confusing.
        if not 'output' in input_no_initial_ns:
            graph.node(f'input_{input}', label=input_no_initial_ns, shape='box', style='rounded, filled', fillcolor='lightgreen')
            graph.edge(f'input_{input}', tool_stem, color=font_edge_color)
    for output in tool_cwl['outputs']:
        output_initial_ns = output.split('___')[0].split('__')[0]
        output_no_initial_ns = output.replace(f'{output_initial_ns}__', '')
        graph.node(f'output_{output}', label=output_no_initial_ns, shape='box', style='rounded, filled', fillcolor='lightyellow')
        graph.edge(tool_stem, f'output_{output}', color=font_edge_color)
    graph.render(format='png')


def make_plugins_DAG(tools: Tools) -> None:
    """Uses the `neato` executable from the graphviz package to make a Directed
    Acyclic Graph consisting of a node for each CWL CommandLineTool and no edges.

    Args:
        tools (Tools): The CWL CommandLineTool definitions found using get_tools_cwl()
    """
    # NOTE: Do not use the default 'dot' engine. Use neato / fdp / sfdp
    # and set pack=0 to remove the massive blank space around each node.
    # Also note that despite my best efforts, I cannot force graphviz to
    # change the aspect ratio (to double the width for making slides)
    # without simply stretching it and distorting the nodes, so we can just
    # partition the tools into two squares and display them side by side.
    num_tools_half = int(len(list(tools)) / 2)
    for i in [0,1]:
        yaml_path = f'autogenerated/DAG/plugins{i}'
        Path(yaml_path).mkdir(parents=True, exist_ok=True)
        graph = graphviz.Digraph(name=yaml_path)
        graph.engine = 'neato'
        graph.attr(pack='0')
        graph.attr(bgcolor="transparent") # Useful for making slides
        graph.attr(fontcolor=font_edge_color)
        for tool in list(tools)[i*num_tools_half:(i+1)*num_tools_half]:
            (tool_path, tool_cwl) = tools[tool]
            graph.node(Path(tool_path).stem, shape='box', style='rounded, filled', fillcolor='lightblue', fontsize="24", width='0.75')
        graph.render(format='png')


def parse_int_string_tuple(string: str) -> Tuple[int, str]:
    # Parses a string of the form '(int, string)'
    string_no_parens = string.strip()[1:-1]
    (str1, str2) = string_no_parens.split(',')
    return (int(str1.strip()), str2.strip())