import argparse
import copy
from pathlib import Path
from typing import List, Tuple

import graphviz
import yaml
from typing import Any, Dict

from . import auto_gen_header
from .wic_types import NodeData, Yaml, RecursiveData, YamlForest

# Use white for dark backgrounds, black for light backgrounds
font_edge_color = 'white'


def step_name_str(yaml_stem: str, i: int, step_key: str) -> str:
    # Use double underscore so we can '__'.split() below.
    # (This should work as long as yaml_stem and step_key do not contain __)
    return f'{yaml_stem}__step__{i+1}__{step_key}'


def parse_step_name_str(step_name_str: str) -> Tuple[str, int, str]:
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


def add_graph_edge(args: argparse.Namespace, graph: graphviz.Digraph, nss1: List[str], nss2: List[str], label: str, color: str = font_edge_color) -> None:
    nss1 = nss1[:(1 + args.graph_inline_depth)]
    edge_node1 = '___'.join(nss1)
    nss2 = nss2[:(1 + args.graph_inline_depth)]
    edge_node2 = '___'.join(nss2)
    # Hide internal self-edges
    if not edge_node1 == edge_node2:
        if args.graph_label_edges:
            graph.edge(edge_node1, edge_node2, color=color, label=label)
        else:
            graph.edge(edge_node1, edge_node2, color=color)


def partition_by_lowest_common_ancestor(nss1: List[str], nss2: List[str]) -> Tuple[List[str], List[str]]:
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
    # Get the dictionary key (i.e. the name) of each step.
    steps_keys = []
    for step in steps:
        steps_keys += list(step)
    #print(steps_keys)
    return steps_keys


def extract_backend(yaml_tree_: Yaml, yaml_path: Path) -> Tuple[str, Yaml]:
    yaml_tree = copy.deepcopy(yaml_tree_)
    backend = ''
    if 'backends' in yaml_tree:
        if 'default_backend' in yaml_tree:
            backend = yaml_tree['default_backend']
            del yaml_tree['default_backend']
        if 'backend' in yaml_tree:
            backend = yaml_tree['backend']
            del yaml_tree['backend']
        if backend == '':
            raise Exception(f'Error! No backend in {yaml_path}!')
        if backend not in yaml_tree['backends']:
            raise Exception(f'Error! No steps for backend {backend} in {yaml_path}!')
        steps = yaml_tree['backends'][backend]['steps']
        del yaml_tree['backends']
        yaml_tree.update({'steps': steps})
    elif 'steps' in yaml_tree:
        pass # steps = yaml_tree['steps']
    else:
        raise Exception(f'Error! No backends and/or steps in {yaml_path}!')
    return (backend, yaml_tree)


def inline_sub_steps(yaml_path: Path, tools: Dict[str, Tuple[str, Any]], yml_paths: Dict[str, Path]) -> List[Yaml]:
    # Load the high-level yaml workflow file.
    with open(Path(yaml_path), 'r') as y:
        yaml_tree: Yaml = yaml.safe_load(y.read())

    (back_name_, yaml_tree) = extract_backend(yaml_tree, yaml_path)
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


def flatten(lists: List[List[Any]]) -> List[List[Any]]:
    return [x for lst in lists for x in lst]


def flatten_recursive_data(recursive_data: RecursiveData) -> List[NodeData]:
    sub_node_data = recursive_data[0]
    sub_rec_lists = [flatten_recursive_data(r) for r in recursive_data[1]]
    return [sub_node_data] + flatten(sub_rec_lists)


def print_forest(forest: YamlForest) -> None:
    print(forest[0][0])
    print(yaml.dump(forest[0][1]))
    print(yaml.dump(forest[1]))


def flatten_forest(forest: YamlForest) -> List[YamlForest]:
    #print_forest(forest)
    if forest == {}:
        return []
    yaml_tree = forest[0][1]
    if 'backends' in yaml_tree:
        #print_forest(forest)
        back_name = ''
        if 'default_backend' in yaml_tree:
            back_name = yaml_tree['default_backend']
        if 'backend' in yaml_tree:
            back_name = yaml_tree['backend']
        if back_name == '':
            print_forest(forest)
            raise Exception(f'Error! No backend in yaml forest!\n')
        step_1 = forest[1][back_name][0][1]['steps'][0]
        step_name_1 = list(step_1.keys())[0]
        if Path(step_name_1).suffix == '.yml':
            # Choose a specific backend
            return flatten_forest(forest[1][back_name])
        else:
            return [forest]

    forests = list(forest[1].values())
    sub_forests = [flatten_forest(f) for f in forests]
    # Use depth first search flattening to match flatten_recursive_data()
    #bfs = forests + flatten(sub_forests)
    dfs_lists = [[f] + fs for f, fs in zip(forests, sub_forests)]
    dfs = flatten(dfs_lists)
    return dfs


def write_to_disk(recursive_data: RecursiveData, path: Path, relative_run_path: bool) -> None:
    sub_node_data = recursive_data[0]
    namespaces = sub_node_data[0]
    yaml_stem = sub_node_data[1]
    yaml_tree = sub_node_data[2]
    yaml_inputs = sub_node_data[3]
    # TODO: destructure sub_node_data rather than using explicit indexing.

    path.mkdir(parents=True, exist_ok=True)
    if relative_run_path:
        filename_cwl = f'{yaml_stem}.cwl'
        filename_yml = f'{yaml_stem}_inputs.yml'
    else:
        filename_cwl = '___'.join(namespaces + [f'{yaml_stem}.cwl'])
        filename_yml = '___'.join(namespaces + [f'{yaml_stem}_inputs.yml'])

    # Dump the compiled yaml file to disk.
    # Use sort_keys=False to preserve the order of the steps.
    yaml_content = yaml.dump(yaml_tree, sort_keys=False, line_break='\n', indent=2)
    with open(path / filename_cwl, 'w') as w:
        w.write('#!/usr/bin/env cwl-runner\n')
        w.write(auto_gen_header)
        w.write(''.join(yaml_content))

    yaml_content = yaml.dump(yaml_inputs, sort_keys=False, line_break='\n', indent=2)
    with open(path / filename_yml, 'w') as inp:
        inp.write(auto_gen_header)
        inp.write(yaml_content)

    for r in recursive_data[1]:
        subpath = path
        if relative_run_path:
            sub_namespace = r[0][0]
            sub_step_name = sub_namespace[-1]
            subpath = path / sub_step_name
        write_to_disk(r, subpath, relative_run_path)


def recursively_delete_dict_key(key: str, obj: Any) -> Any:
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
            