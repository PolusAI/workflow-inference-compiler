import argparse
from pathlib import Path
from typing import Dict, List

from . import utils
from .wic_types import (GraphReps, InternalOutputs, Namespaces, StepId, Tools,
                        WorkflowOutputs, Yaml)


def maybe_add_subworkflow_requirement(yaml_tree: Yaml, tools: Tools, steps_keys: List[str],
                                      wic_steps: Yaml, subkeys: List[str]) -> None:
    """Adds a SubworkflowFeatureRequirement if there are any subworkflows

    Args:
        yaml_tree (Yaml): A tuple of name and yml AST
        tools (Tools): The CWL CommandLineTool definitions found using get_tools_cwl()
        steps_keys (List[str]): The name of each step in the current CWL workflow
        wic_steps (Yaml): The metadata assocated with the workflow steps
        subkeys (List[str]): The keys associated with subworkflows
    """
    # If there is at least one subworkflow, add a SubworkflowFeatureRequirement
    subs = []
    for i, step_key in enumerate(steps_keys):
        sub_wic = wic_steps.get(f'({i+1}, {step_key})', {})
        plugin_ns_i = sub_wic.get('wic', {}).get('namespace', 'global')
        if step_key not in subkeys:
            step_id = StepId(Path(step_key).stem, plugin_ns_i)
            sub = tools[step_id].cwl['class'] == 'Workflow'
            subs.append(sub)
    if (not subkeys == []) or any(subs):
        subworkreq = 'SubworkflowFeatureRequirement'
        subworkreqdict = {subworkreq: {'class': subworkreq}}
        if 'requirements' in yaml_tree:
            if not subworkreq in yaml_tree['requirements']:
                new_reqs = dict(list(yaml_tree['requirements'].items()) + list(subworkreqdict))
                yaml_tree['requirements'].update(new_reqs)
        else:
            yaml_tree['requirements'] = subworkreqdict


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
    """add_yamldict_keyval partially applied with in_out='in' """
    return add_yamldict_keyval(steps_i, step_key, 'in', keyval)

def add_yamldict_keyval_out(steps_i: Yaml, step_key: str, keyval: List[str]) -> Yaml:
    """add_yamldict_keyval partially applied with in_out='out' """
    return add_yamldict_keyval(steps_i, step_key, 'out', keyval) # type: ignore


def get_workflow_outputs(args: argparse.Namespace,
                         namespaces: Namespaces,
                         is_root: bool,
                         yaml_stem: str,
                         steps: List[Yaml],
                         wic_steps: Yaml,
                         outputs_workflow: WorkflowOutputs,
                         vars_workflow_output_internal: InternalOutputs,
                         graph: GraphReps,
                         tools: Tools,
                         step_node_name: str) -> Dict[str, Dict[str, str]]:
    """Chooses a subset of the CWL outputs: to actually output

    Args:
        args (argparse.Namespace): The command line arguments
        namespaces (Namespaces): Specifies the path in the AST of the current subworkflow
        is_root (bool): True if this is the root workflow
        yaml_stem (str): The name of the current subworkflow (stem of the yaml filepath)
        steps (List[Yaml]): The steps: tag of a CWL workflow
        wic_steps (Yaml): The metadata assocated with the workflow steps
        outputs_workflow (WorkflowOutputs): Contains the contents of the out: tags for each step.
        vars_workflow_output_internal (InternalOutputs): Keeps track of output\n
        variables which are internal to the root workflow, but not necessarily to subworkflows.
        graph (GraphReps): A tuple of a GraphViz DiGraph and a networkx DiGraph
        tools (Tools): The CWL CommandLineTool definitions found using get_tools_cwl()
        step_node_name (str): The namespaced name of the current step

    Returns:
        Dict[str, Dict[str, str]]: The actual outputs to be specified in the generated CWL file
    """
    # Add the outputs of each step to the workflow outputs
    workflow_outputs = {}
    steps_keys = utils.get_steps_keys(steps)
    for i, step_key in enumerate(steps_keys):
        sub_wic = wic_steps.get(f'({i+1}, {step_key})', {})
        plugin_ns_i = sub_wic.get('wic', {}).get('namespace', 'global')
        stem = Path(step_key).stem
        step_id = StepId(stem, plugin_ns_i)
        tool_i = tools[step_id].cwl
        step_name_i = utils.step_name_str(yaml_stem, i, step_key)
        out_keys = steps[i][step_key]['out']
        for out_key in out_keys:
            out_var = f'{step_name_i}/{out_key}'
            # Avoid duplicating intermediate outputs in GraphViz
            out_key_no_namespace = out_key.split('___')[-1]
            if args.graph_show_outputs:
                vars_nss = [var.replace('/', '___') for var in vars_workflow_output_internal]
                case1 = (tool_i['class'] == 'Workflow') and (not out_key in vars_nss)
                # Avoid duplicating outputs from subgraphs in parent graphs.
                output_node_name = '___'.join(namespaces + [step_name_i, out_key])
                lengths_off_by_one = (len(step_node_name.split('___')) + 1 == len(output_node_name.split('___')))
                # TODO: check is_root here
                case1 = case1 and not is_root and lengths_off_by_one
                case2 = (tool_i['class'] == 'CommandLineTool') and (not out_var in vars_workflow_output_internal)
                if case1 or case2:
                    graph_gv = graph.graphviz
                    graph_nx = graph.networkx
                    graphdata = graph.graphdata
                    attrs = {'label': out_key_no_namespace, 'shape': 'box',
                             'style': 'rounded, filled', 'fillcolor': 'lightyellow'}
                    graph_gv.node(output_node_name, **attrs)
                    font_edge_color = 'black' if args.graph_dark_theme else 'white'
                    if args.graph_label_edges:
                        graph_gv.edge(step_node_name, output_node_name, color=font_edge_color,
                                      label=out_key_no_namespace)  # Is labeling necessary?
                    else:
                        graph_gv.edge(step_node_name, output_node_name, color=font_edge_color)
                    graph_nx.add_node(output_node_name)
                    graph_nx.add_edge(step_node_name, output_node_name)
                    graphdata.nodes.append((output_node_name, attrs))
                    graphdata.edges.append((step_node_name, output_node_name, {}))
            # NOTE: Unless we are in the root workflow, we always need to
            # output everything. This is because while we are within a
            # subworkflow, we do not yet know if a subworkflow output will be used as
            # an input in a parent workflow (either explicitly, or using inference).
            # Thus, we have to output everything.
            # However, once we reach the root workflow, we can choose whether
            # we want to pollute our output directory with intermediate files.
            # (You can always enable --provenance to get intermediate files.)
            # NOTE: Remove is_root for now because in test_cwl_embedding_independence,
            # we recompile all subworkflows as if they were root. Thus, for now
            # we need to enable args.cwl_output_intermediate_files
            # Exclude intermediate 'output' files.
            if out_var in vars_workflow_output_internal and not args.cwl_output_intermediate_files: # and is_root
                continue
            out_name = f'{step_name_i}___{out_key}'  # Use triple underscore for namespacing so we can split later
            #print('out_name', out_name)

        for out_key, out_dict in outputs_workflow[i].items():
            out_name = f'{step_name_i}___{out_key}'  # Use triple underscore for namespacing so we can split later
            out_var = f'{step_name_i}/{out_key}'
            workflow_outputs.update({out_name: {**out_dict, 'outputSource': out_var}})
        #print('workflow_outputs', workflow_outputs)
    return workflow_outputs
