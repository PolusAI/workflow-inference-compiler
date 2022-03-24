from pathlib import Path
import yaml

# Use white for dark backgrounds, black for light backgrounds
font_edge_color = 'white'


def add_yamldict_keyval(steps_i, step_key, in_out, keyval):
    if steps_i[step_key]:
        if in_out in steps_i[step_key]:
            new_keys = dict(list(steps_i[step_key][in_out].items()) + list(keyval.items()))
            new_keyvals = [(k, v) if k != in_out else (k, new_keys) for k, v in steps_i[step_key].items()]
        else:
            new_keys = keyval
            new_keyvals = dict(list(steps_i[step_key].items()) + [(in_out, new_keys)])
        steps_i[step_key].update(new_keyvals)
    else:
        steps_i = {step_key: {in_out: keyval}}
    return steps_i


def add_graph_edge(args, graph, nss1, nss2, label):
    nss1 = nss1[:(1 + args.graph_inline_depth)]
    edge_node1 = '___'.join(nss1)
    nss2 = nss2[:(1 + args.graph_inline_depth)]
    edge_node2 = '___'.join(nss2)
    # Hide internal self-edges
    if not edge_node1 == edge_node2:
        if args.graph_label_edges:
            graph.edge(edge_node1, edge_node2, color=font_edge_color, label=label)
        else:
            graph.edge(edge_node1, edge_node2, color=font_edge_color)


def partition_by_lowest_common_ancestor(nss1, nss2):
    # Only partition nss1; if you want to partition nss1
    # just switch the arguments at the call site.
    if nss1 == [] or nss2 == []:
        return ([], nss1) # Base case
    if nss1[0] == nss2[0]: # Keep going
        (nss1_heads, nss1_tails) = partition_by_lowest_common_ancestor(nss1[1:], nss2[1:])
        return ([nss1[0]] + nss1_heads, nss1_tails)
    else:
        return ([], nss1)


def extract_backend_steps(yaml_tree, yaml_path):
    backend = None
    if 'backends' in yaml_tree:
        if 'default_backend' in yaml_tree:
            backend = yaml_tree['default_backend']
            del yaml_tree['default_backend']
        if backend is None:
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
    return yaml_tree


def inline_sub_steps(yaml_path, tools, yml_paths):
    # Load the high-level yaml workflow file.
    with open(Path(yaml_path), 'r') as y:
        yaml_tree = yaml.safe_load(y.read())

    yaml_tree = extract_backend_steps(yaml_tree, yaml_path)
    steps = yaml_tree['steps']

    # Get the dictionary key (i.e. the name) of each step.
    steps_keys = []
    for step in steps:
        steps_keys += list(step)

    subkeys = [key for key in steps_keys if key not in tools]

    steps_all = []
    for i, step_key in enumerate(steps_keys):
        if step_key in subkeys:
            path = Path(yml_paths[Path(step_key).stem])
            steps_i = inline_sub_steps(path, tools, yml_paths)
        else:
            steps_i = [steps[i]]
        steps_all.append(steps_i)

    steps_all_flattened = [step for steps in steps_all for step in steps]
    return steps_all_flattened