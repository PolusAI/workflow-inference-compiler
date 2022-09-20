import argparse
import copy
from pathlib import Path
from typing import Any, Dict, List

from . import utils
from .wic_types import (GraphReps, InternalOutputs, Namespaces, StepId, Tool, Tools,
                        WorkflowOutputs, Yaml)


def maybe_add_requirements(yaml_tree: Yaml, tools: Tools, steps_keys: List[str],
                           wic_steps: Yaml, subkeys: List[str]) -> None:
    """Adds any necessary CWL requirements

    Args:
        yaml_tree (Yaml): A tuple of name and yml AST
        tools (Tools): The CWL CommandLineTool definitions found using get_tools_cwl()
        steps_keys (List[str]): The name of each step in the current CWL workflow
        wic_steps (Yaml): The metadata assocated with the workflow steps
        subkeys (List[str]): The keys associated with subworkflows
    """
    bools = []
    subwork = []
    scatter = []
    stepinp = []
    for i, step_key in enumerate(steps_keys):
        sub_wic = wic_steps.get(f'({i+1}, {step_key})', {})

        if step_key not in subkeys:
            plugin_ns_i = sub_wic.get('wic', {}).get('namespace', 'global')
            step_id = StepId(Path(step_key).stem, plugin_ns_i)
            sub = tools[step_id].cwl['class'] == 'Workflow'
            bools.append(sub)

        if 'scatter' in sub_wic:
            scatter = ['ScatterFeatureRequirement']

        in_step = yaml_tree['steps'][i][step_key].get('in')
        sub_wic_copy = copy.deepcopy(sub_wic)
        if 'wic' in sub_wic_copy:
            del sub_wic_copy['wic']
        if (utils.recursively_contains_dict_key('valueFrom', in_step) or
            utils.recursively_contains_dict_key('valueFrom', sub_wic_copy)):
            stepinp = ['StepInputExpressionRequirement', 'InlineJavascriptRequirement']

    if (not subkeys == []) or any(bools):
        subwork = ['SubworkflowFeatureRequirement']

    reqs = subwork + scatter + stepinp
    if reqs:
        reqsdict: Dict[str, Dict] = {r: {} for r in reqs}
        if 'requirements' in yaml_tree:
            new_reqs = dict(list(yaml_tree['requirements'].items()) + list(reqsdict))
            yaml_tree['requirements'].update(new_reqs)
        else:
            yaml_tree['requirements'] = reqsdict


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
                         tools_lst: List[Tool],
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
        tools_lst (List[Tool]): A list of the CWL CommandLineTools or compiled subworkflows for the current workflow.
        step_node_name (str): The namespaced name of the current step

    Returns:
        Dict[str, Dict[str, str]]: The actual outputs to be specified in the generated CWL file
    """
    # Add the outputs of each step to the workflow outputs
    workflow_outputs = {}
    steps_keys = utils.get_steps_keys(steps)
    for i, step_key in enumerate(steps_keys):
        tool_i = tools_lst[i].cwl
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

        wic_step_i = wic_steps.get(f'({i+1}, {step_key})', {})
        for out_key, out_dict in outputs_workflow[i].items():
            if 'scatter' in wic_step_i:
                out_dict['type'] = canonicalize_type(out_dict['type'])
                # Promote scattered output types to arrays
                out_dict['type'] = {'type': 'array', 'items': out_dict['type']}

            out_name = f'{step_name_i}___{out_key}'  # Use triple underscore for namespacing so we can split later
            out_var = f'{step_name_i}/{out_key}'
            workflow_outputs.update({out_name: {**out_dict, 'outputSource': out_var}})
        #print('workflow_outputs', workflow_outputs)
    # NOTE: The fix_conflicts 'feature' of cwltool prevents files from being
    # overwritten by appending _2, _3 etc.
    # The problem is that these renamed files now no longer match the glob
    # patterns in the outputBinding tags, thus they are not copied to the
    # final output folder in relocateOutputs() and/or stage_files().
    # Note that this error is not detected using --validate.
    # One workaround is to simply output all files.
    # TODO: glob "." is still returning null; need to use InitialWorkDirRequirement??
    output_all = {'output_all':
                    {'type':
                        {'type': 'array',
                         'items': ['Directory', 'File']},
                     'outputBinding': {'glob': '\".\"'},
                     'format': 'edam:format_2330'}} # 'Textual format'
    workflow_outputs.update(output_all) # type: ignore
    return workflow_outputs


def canonicalize_type(type_obj: Any) -> Any:
    """Recursively desugars the CWL type: field into a canonical normal form.\n
    In particular, CWL automatically desugars File[] into {'type': 'array', 'items': File},
    but File[][] causes a syntax error! Etc.

    Args:
        type_obj (Any): An object that is a syntactic hodgepodge of valid CWL types.

    Returns:
        Any: The JSON canonical normal form associated with type_obj
    """
    if isinstance(type_obj, str):
        if len(type_obj) >= 1 and type_obj[-1:] == '?':
            return ['null', canonicalize_type(type_obj[:-1])]
        if len(type_obj) >= 2 and type_obj[-2:] == '[]':
            return {'type': 'array', 'items': canonicalize_type(type_obj[:-2])}
    if isinstance(type_obj, Dict):
        if type_obj.get('type') == 'array':
            return {**type_obj, 'items': canonicalize_type(type_obj['items'])}
    return type_obj
            