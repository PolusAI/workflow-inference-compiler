import argparse
import copy
from pathlib import Path
import subprocess as sub
from typing import Dict, List

import graphviz
import networkx  as nx

from . import inference
from . import utils
from .wic_types import Yaml, Tools, DollarDefs, DollarCalls, CompilerInfo, WorkflowInputsFile, InternalOutputs, Graph, YamlTree

# Use white for dark backgrounds, black for light backgrounds
font_edge_color = 'white'


def compile_workflow(yaml_tree_: YamlTree,
                     args: argparse.Namespace,
                     namespaces: List[str],
                     subgraphs: List[Graph],
                     vars_dollar_defs: DollarDefs,
                     vars_dollar_calls: DollarCalls,
                     tools: Tools,
                     is_root: bool,
                     relative_run_path: bool) -> CompilerInfo:
    (yaml_path, yaml_tree) = yaml_tree_

    # Check for top-level yml dsl args
    wic = {'wic': yaml_tree.get('wic', {})}
    #print('wic', yaml.dump(wic))
    if 'wic' in yaml_tree:
        del yaml_tree['wic']
    wic_steps = wic['wic'].get('steps', {})
    
    yaml_stem = Path(yaml_path).stem

    (back_name_, yaml_tree) = utils.extract_backend(yaml_tree, Path(yaml_path))
    steps: List[Yaml] = yaml_tree['steps']

    steps_keys = utils.get_steps_keys(steps)

    subkeys = [key for key in steps_keys if key not in tools]

    # Add headers
    yaml_tree['cwlVersion'] = 'v1.0'
    yaml_tree['class'] = 'Workflow'
    yaml_tree['$namespaces'] = {'edam': 'https://edamontology.org/'}
    yaml_tree['$schemas'] = ['https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl']

    # NOTE: currently mutates yaml_tree (maybe)
    maybe_add_subworkflow_requirement(yaml_tree, tools, steps_keys, subkeys)

    # Collect workflow input parameters
    inputs_workflow = {}
    inputs_file_workflow = {}
    
    # Collect the internal workflow output variables
    vars_workflow_output_internal = []

    # Collect recursive dollar variable definitions.
    vars_dollar_defs_copy = copy.deepcopy(vars_dollar_defs)

    # Collect recursive dollar variable call sites.
    vars_dollar_calls_copy = copy.deepcopy(vars_dollar_calls)

    # Collect recursive subworkflow data
    step_1_names = []
    sibling_subgraphs = []
    
    recursive_data_list = []

    graph = subgraphs[-1] # Get the current graph
    (graph_gv, graph_nx) = graph

    # TODO: Check for top-level yml dsl args

    for i, step_key in enumerate(steps_keys):
        step_name_i = utils.step_name_str(yaml_stem, i, step_key)
        stem = Path(step_key).stem

        # Recursively compile subworkflows, adding compiled cwl file contents to tools
        if step_key in subkeys:
            # Extract the sub yaml file that we pre-loaded from disk.
            sub_yaml_tree = (step_key, steps[i][step_key])

            subgraph_gv = graphviz.Digraph(name=f'cluster_{step_key}')
            subgraph_gv.attr(label=step_key) # str(path)
            subgraph_gv.attr(color='lightblue')  # color of outline
            subgraph_nx = nx.DiGraph()
            subgraph = (subgraph_gv, subgraph_nx)

            steps[i].update({step_key: {}}) # delete yml subtree
            subworkflow_data = compile_workflow(sub_yaml_tree, args, namespaces + [step_name_i], subgraphs + [subgraph], vars_dollar_defs, vars_dollar_calls, tools, False, relative_run_path)
            
            recursive_data = subworkflow_data[0]
            recursive_data_list.append(recursive_data)

            sub_node_data = recursive_data[0]
            # TODO: destructure sub_node_data rather than using explicit indexing.
            
            sibling_subgraphs.append(sub_node_data[-1]) # TODO: Just subgraph?
            step_1_names.append(subworkflow_data[-1])
            # Add compiled cwl file contents to tools
            tools[stem] = (stem + '.cwl', sub_node_data[2])

            # Initialize the above from recursive values.
            # Do not initialize inputs_workflow. See comment below.
            # inputs_workflow.update(subworkflow_data[1])
            inputs_namespaced = dict([(f'{step_name_i}___{k}', val) for k, val in subworkflow_data[2].items()]) # _{step_key}_input___{k}
            inputs_file_workflow.update(inputs_namespaced)
            vars_workflow_output_internal += subworkflow_data[3]
            vars_dollar_defs.update(subworkflow_data[4])
            vars_dollar_calls.update(subworkflow_data[5])

        tool_i = tools[stem][1]
        utils.make_tool_DAG(stem, tools[stem])

        # Add run tag, using relative or flat-directory paths
        # NOTE: run: path issues were causing test_cwl_embedding_independence()
        # to fail, so I simply ignore the run tag in that test.
        run_path = tools[stem][0]
        if relative_run_path:
            if step_key in subkeys:
                run_path = step_name_i + '/' + run_path
            else:
                run_path = ('../' * len(namespaces)) + run_path
                run_path = '../' + run_path # Relative to autogenerated/
        else:
            if step_key in subkeys:
                run_path = '___'.join(namespaces + [step_name_i, run_path])
            else:
                run_path = '../' + run_path # Relative to autogenerated/

        if steps[i][step_key]:
            if not 'run' in steps[i][step_key]:
                steps[i][step_key].update({'run': run_path})
        else:
            steps[i] = {step_key: {'run': run_path}}

        # Generate intermediate file names between steps.
        # The inputs for a given step need to come from either
        # 1. the input.yaml file for the overall workflow (extract into separate yml file) or
        # 2. the outputs from the previous steps (most recent first).
        # If there isn't an exact match, remove input_* and output_* from the
        # current step and previous steps, respectively, and then check again.
        # If there still isn't an exact match, explicit renaming may be required.

        args_provided = []
        if steps[i][step_key] and 'in' in steps[i][step_key]:
            args_provided = list(steps[i][step_key]['in'])
        #print(args_provided)

        in_tool = tool_i['inputs']
        #print(list(in_tool.keys()))
        if tool_i['class'] == 'CommandLineTool':
            args_required = [arg for arg in in_tool if not (in_tool[arg].get('default') or in_tool[arg]['type'][-1] == '?')]
        elif tool_i['class'] == 'Workflow':
            args_required = [arg for arg in in_tool]
            
            # Add the inputs. For now, assume that all sub-workflows have been
            # auto-generated from a previous application of this script, and
            # thus that all of their transitive inputs have been satisfied.
            # (i.e. simply combine the input yml files using the cat command.)
            steps[i][step_key]['in'] = dict([(key, key) for key in args_required])
        else:
            raise Exception(f'Unknown class', tool_i['class'])
        
        # Note: Some config tags are not required in the cwl files, but are in
        # fact required in the python source code! See check_mandatory_property
        # (Solution: refactor all required arguments out of config and list
        # them as explicit inputs in the cwl files, then modify the python
        # files accordingly.)
        #print(args_required)

        sub_args_provided = [arg for arg in args_required if arg in vars_dollar_calls]
        #print(sub_args_provided)

        label = step_key
        if args.graph_label_stepname:
            label = step_name_i
        step_node_name = '___'.join(namespaces + [step_name_i])
        if not tool_i['class'] == 'Workflow':
            graph_gv.node(step_node_name, label=label, shape='box', style='rounded, filled', fillcolor='lightblue')
            graph_nx.add_node(step_node_name)
        elif not (step_key in subkeys and len(namespaces) < args.graph_inline_depth):
            nssnode = namespaces + [step_name_i]
            nssnode = nssnode[:(1 + args.graph_inline_depth)]
            step_node_name = '___'.join(nssnode)
            graph_gv.node(step_node_name, label=label, shape='box', style='rounded, filled', fillcolor='lightblue')
            graph_nx.add_node(step_node_name)

        # NOTE: sub_args_provided are handled within the args_required loop below
        for arg_key in args_provided:
            # Extract input value into separate yml file
            # Replace it here with a new variable name
            arg_val = steps[i][step_key]['in'][arg_key]
            in_name = f'{step_name_i}___{arg_key}'  # Use triple underscore for namespacing so we can split later # {step_name_i}_input___{arg_key}
            in_type = in_tool[arg_key]['type'].replace('?', '')  # Providing optional arguments makes them required
            in_dict = {'type': in_type}
            if 'format' in in_tool[arg_key]:
                in_format = in_tool[arg_key]['format']
                in_dict['format'] = in_format
            if arg_val[0] == '&':
                arg_val = arg_val[1:]  # Remove &
                #print('arg_key, arg_val', arg_key, arg_val)
                # NOTE: There can only be one definition, but multiple call sites.
                if not vars_dollar_defs.get(arg_val):
                    # if first time encountering arg_val, i.e. if defining
                    inputs_workflow.update({in_name: in_dict})
                    in_dict = {**in_dict, 'value': arg_val}
                    inputs_file_workflow.update({in_name: in_dict})
                    steps[i][step_key]['in'][arg_key] = in_name
                    vars_dollar_defs.update({arg_val: (namespaces + [step_name_i], arg_key)})
                    # TODO: Show input node?
                else:
                    raise Exception(f"Error! Multiple definitions of &{arg_val}!")
            elif arg_val[0] == '*':
                arg_val = arg_val[1:]  # Remove *
                if not vars_dollar_defs.get(arg_val):
                    if is_root:
                        # TODO: Check this comment.
                        # Even if is_root, we don't want to raise an Exception
                        # here because in test_cwl_embedding_independence, we
                        # recompile all subworkflows as if they were root. That
                        # will cause this code path to be taken but it is not
                        # actually an error. Add a CWL input for testing only.
                        print(f"Error! No definition found for &{arg_val}!")
                        print(f"Creating the CWL input {in_name} anyway, but")
                        print("without any corresponding input value this will fail validation!")
                    inputs_workflow.update({in_name: in_dict})
                    steps[i][step_key]['in'][arg_key] = in_name
                else:
                    (nss_def_init, var) =  vars_dollar_defs[arg_val]

                    nss_def_embedded = var.split('___')[:-1]
                    nss_call_embedded = arg_key.split('___')[:-1]
                    nss_def = nss_def_init + nss_def_embedded
                    # [step_name_i] is correct; nss_def_init already contains step_name_j from the recursive call
                    nss_call = namespaces + [step_name_i] + nss_call_embedded

                    nss_def_inits, nss_def_tails = utils.partition_by_lowest_common_ancestor(nss_def, nss_call)
                    nss_call_inits, nss_call_tails = utils.partition_by_lowest_common_ancestor(nss_call, nss_def)
                    # nss_def and nss_call are paths into the abstract 'call stack'.
                    # This defines the 'common node' in the call stack w.r.t. the inits.
                    assert nss_def_inits == nss_call_inits
                    
                    # TODO: Check this comment.
                    # Relative to the common node, if the call site of an explicit
                    # edge is at a depth > 1, (i.e. if it is NOT simply of the form
                    # last_namespace/input_variable) then we
                    # need to create inputs in all of the intervening CWL files
                    # so we can pass in the values from the outer scope(s). Here,
                    # we simply need to use in_name and add to inputs_workflow
                    # and vars_dollar_calls. The outer scope(s) are handled by
                    # the sub_args_provided clause below.
                    # Note that the outputs from the definition site are bubbled
                    # up the call stack until they reach the common node.
                    if len(nss_call_tails) > 1:
                        inputs_workflow.update({in_name: in_dict})
                        steps[i][step_key]['in'][arg_key] = in_name
                        # Store var_dollar call site info up through the recursion.
                        vars_dollar_calls.update({in_name: vars_dollar_defs[arg_val]}) # {in_name, (namespaces + [step_name_i], var)} ?
                    elif len(nss_call_tails) == 1:
                        # TODO: Check this comment.
                        # The definition site recursion (only, if any) has completed
                        # and we are already in the common node, thus
                        # we need to pass in the value from the definition site.
                        # Note that since len(nss_call_tails) == 1,
                        # there will not be any call site recursion in this case.
                        var_slash = nss_def_tails[0] + '/' + '___'.join(nss_def_tails[1:] + [var])
                        steps[i][step_key]['in'][arg_key] = var_slash
                    else:
                        # Since nss_call includes step_name_i, this should never happen...
                        raise Exception("Error! len(nss_call_tails) == 0! Please file a bug report!\n" +
                                        f'nss_def {nss_def}\n nss_call {nss_call}')

                    # Add an edge, but in a carefully chosen subgraph.
                    # If you add an edge whose head/tail is outside of the subgraph,
                    # graphviz may segfault! Moreover, even if graphviz doesn't
                    # segfault, adding an edge in a given subgraph can cause the
                    # nodes themselves to be rendered in that subgraph, even
                    # though the nodes are defined in a different subgraph!
                    # The correct thing to do is to use the graph associated with
                    # the lowest_common_ancestor of the definition and call site.
                    # (This is the only reason we need to pass in all subgraphs.)
                    label = var.split('___')[-1]
                    graph_init = subgraphs[len(nss_def_inits)]
                    # Let's use regular blue for explicit edges.
                    # Use constraint=false ?
                    utils.add_graph_edge(args, graph_init, nss_def, nss_call, label, color='blue')
            else:
                inputs_workflow.update({in_name: in_dict})
                in_dict = {**in_dict, 'value': arg_val}
                inputs_file_workflow.update({in_name: in_dict})
                steps[i][step_key]['in'][arg_key] = in_name
                if args.graph_show_inputs:
                    input_node_name = '___'.join(namespaces + [step_name_i, arg_key])
                    graph_gv.node(input_node_name, label=arg_key, shape='box', style='rounded, filled', fillcolor='lightgreen')
                    graph_gv.edge(input_node_name, step_node_name, color=font_edge_color)
                    graph_nx.add_node(input_node_name)
                    graph_nx.add_edge(input_node_name, step_node_name)
        
        for arg_key in args_required:
            #print('arg_key', arg_key)
            in_name = f'{step_name_i}___{arg_key}'
            if arg_key in args_provided:
                continue  # We already covered this case above.
            if arg_key in sub_args_provided: # Edges have been explicitly provided
                # The definition site recursion (if any) and the call site
                # recursion (yes, see above), have both completed and we are
                # now in the common node, thus 
                # we need to pass in the value from the definition site.
                # Extract the stored defs namespaces from vars_dollar_calls.
                # (See massive comment above.)
                (nss_def_init, var) = vars_dollar_calls[arg_key]

                nss_def_embedded = var.split('___')[:-1]
                nss_call_embedded = arg_key.split('___')[:-1]
                nss_def = nss_def_init + nss_def_embedded
                # [step_name_i] is correct; nss_def_init already contains step_name_j from the recursive call
                nss_call = namespaces + [step_name_i] + nss_call_embedded

                nss_def_inits, nss_def_tails = utils.partition_by_lowest_common_ancestor(nss_def, nss_call)
                nss_call_inits, nss_call_tails = utils.partition_by_lowest_common_ancestor(nss_call, nss_def)
                assert nss_def_inits == nss_call_inits

                nss_call_tails_stems = [utils.parse_step_name_str(x)[0] for x in nss_call_tails]
                if yaml_stem in nss_call_tails_stems and nss_call_tails_stems.index(yaml_stem) > 0:
                    # i.e. if it is possible to do more recursion
                    # NOTE: This works, and test_cwl_embedding_independence()
                    # passes, but it is NOT morally embedding independent.
                    in_type = in_tool[arg_key]['type']
                    in_dict = {'type': in_type}
                    if 'format' in in_tool[arg_key]:
                        in_format = in_tool[arg_key]['format']
                        in_dict['format'] = in_format
                    inputs_workflow.update({in_name: in_dict})
                    steps[i][step_key]['in'][arg_key] = in_name
                    # Store var_dollar call site info up through the recursion.
                    vars_dollar_calls.update({in_name: vars_dollar_calls[arg_key]})
                else:
                    # TODO: Check this comment.
                    # The definition site recursion (only, if any) has completed
                    # and we are already in the common node, thus
                    # we need to pass in the value from the definition site.
                    # Note that since len(nss_call_tails) == 1,
                    # there will not be any call site recursion in this case.
                    var_slash = nss_def_tails[0] + '/' + '___'.join(nss_def_tails[1:] + [var])
                    steps[i][step_key]['in'][arg_key] = var_slash

                # NOTE: We already added an edge to the appropriate subgraph above.
                # TODO: vars_workflow_output_internal?
            else:
                in_name_in_inputs_file_workflow: bool = (in_name in inputs_file_workflow)
                steps[i] = inference.perform_edge_inference(args, tools, steps_keys,
                    yaml_stem, i, steps[i], arg_key, graph, is_root, namespaces,
                    vars_workflow_output_internal, inputs_workflow, in_name_in_inputs_file_workflow)
                # NOTE: For now, perform_edge_inference mutably appends to
                # inputs_workflow and vars_workflow_output_internal.
        
        # Add CommandLineTool outputs tags to workflow out tags.
        # Note: Add all output tags for now, but depending on config options,
        # not all output files will be generated. This may cause an error.
        out_keys = list(tool_i['outputs'])
        #print('out_keys', out_keys)
        steps[i] = utils.add_yamldict_keyval_out(steps[i], step_key, out_keys)

        #print()

    # NOTE: add_subgraphs currently mutates graph
    gv_options = wic['wic'].get('graphviz', {})
    ranksame_strs = gv_options.get('ranksame', [])
    import ast
    ranksame_pairs = [ast.literal_eval(x) for x in ranksame_strs]
    steps_ranksame = ['___'.join(namespaces + [utils.step_name_str(yaml_stem, num-1, name)]) for num, name in ranksame_pairs]
    steps_ranksame = [f'"{x}"' for x in steps_ranksame]  # Escape with double quotes.
    add_subgraphs(args, graph, sibling_subgraphs, namespaces, step_1_names, steps_ranksame)
    step_name_1 = get_step_name_1(step_1_names, yaml_stem, namespaces, steps_keys, subkeys)

    # Add the provided inputs of each step to the workflow inputs
    yaml_tree.update({'inputs': inputs_workflow})

    vars_workflow_output_internal = list(set(vars_workflow_output_internal))  # Get uniques
    # (Why are we getting uniques?)
    workflow_outputs = get_workflow_outputs(args, namespaces, is_root, yaml_stem, steps, vars_workflow_output_internal, graph, tools, step_node_name)
    yaml_tree.update({'outputs': workflow_outputs})

    # Finally, rename the steps to be unique
    # and convert the list of steps into a dict
    steps_dict = {}
    for i, step_key in enumerate(steps_keys):
        step_name_i = utils.step_name_str(yaml_stem, i, step_key)
        #steps[i] = {step_name_i: steps[i][step_key]}
        steps_dict.update({step_name_i: steps[i][step_key]})
    yaml_tree.update({'steps': steps_dict})

    # Dump the workflow inputs to a separate yml file.
    yaml_inputs: WorkflowInputsFile = {}
    for key, in_dict in inputs_file_workflow.items():
        new_keyval: WorkflowInputsFile = {}
        if 'File' == in_dict['type']:
            in_format = in_dict['format']
            if isinstance(in_format, List):
                print(f'NOTE: More than one input file format for {key}')
                print(f'formats: {in_format}')
                print(f'Choosing {in_format[0]}')
                in_format = in_format[0]
            new_keyval = {key: {'class': 'File', 'path': in_dict['value'], 'format': in_format}}
        elif 'string' == in_dict['type']:
            # We cannot store string values as a dict, so use type: ignore
            new_keyval = {key: in_dict['value']} # type: ignore 
        else:
            raise Exception(f"Error! Unknown type: {in_dict['type']}")
        yaml_inputs.update(new_keyval)

    if args.cwl_validate:
        print(f'Validating {yaml_stem}.cwl ...')
        cmd = ['cwltool', '--validate', f'{yaml_stem}.cwl']
        sub.run(cmd)
    
    # Note: We do not necessarily need to return inputs_workflow.
    # 'Internal' inputs are encoded in yaml_tree. See Comment above.
    node_data = (namespaces, yaml_stem, yaml_tree, yaml_inputs, vars_dollar_defs_copy, vars_dollar_calls_copy, graph) # step_name_1, plugin_id?
    return ((node_data, recursive_data_list), inputs_workflow, inputs_file_workflow, vars_workflow_output_internal, vars_dollar_defs, vars_dollar_calls, step_name_1)


def maybe_add_subworkflow_requirement(yaml_tree: Yaml, tools: Tools, steps_keys: List[str], subkeys: List[str]) -> None:
    # If there is at least one subworkflow, add a SubworkflowFeatureRequirement
    if (not subkeys == []) or any([tools[Path(key).stem][1]['class'] == 'Workflow' for key in steps_keys if key not in subkeys]):
        subworkreq = 'SubworkflowFeatureRequirement'
        subworkreqdict = {subworkreq: {'class': subworkreq}}
        if 'requirements' in yaml_tree:
            if not subworkreq in yaml_tree['requirements']:
                new_reqs = dict(list(yaml_tree['requirements'].items()) + list(subworkreqdict))
                yaml_tree['requirements'].update(new_reqs)
        else:
            yaml_tree['requirements'] = subworkreqdict


def add_subgraphs(args: argparse.Namespace,
                  graph: Graph,
                  sibling_subgraphs: List[Graph],
                  namespaces: List[str],
                  step_1_names: List[str],
                  steps_ranksame: List[str]) -> None:
    (graph_gv, graph_nx) = graph
    # Add the cluster subgraphs to the main graph, but we need to add them in
    # reverse order to trick the graphviz layout algorithm.
    for sibling in sibling_subgraphs[::-1]: # Reverse!
        (sib_graph_gv, sib_graph_nx) = sibling
        if len(namespaces) < args.graph_inline_depth:
            graph_gv.subgraph(sib_graph_gv)
        graph_nx.add_nodes_from(sib_graph_nx.nodes)
        graph_nx.add_edges_from(sib_graph_nx.edges)
    # Align the cluster subgraphs using the same rank as the first node of each subgraph.
    # See https://stackoverflow.com/questions/6824431/placing-clusters-on-the-same-rank-in-graphviz
    if len(namespaces) < args.graph_inline_depth:
        step_1_names_display = [name for name in step_1_names if len(name.split('___')) < 2 + args.graph_inline_depth]
        if len(step_1_names_display) > 1:
            nodes_same_rank = '\t{rank=same; ' + '; '.join(step_1_names_display) + '}\n'
            graph_gv.body.append(nodes_same_rank)
        if len(steps_ranksame) > 1:
            nodes_same_rank = '\t{rank=same; ' + '; '.join(steps_ranksame) + '}\n'
            graph_gv.body.append(nodes_same_rank)
        


def get_step_name_1(step_1_names: List[str],
                    yaml_stem: str,
                    namespaces: List[str],
                    steps_keys: List[str],
                    subkeys: List[str]) -> str:
    if steps_keys[0] in subkeys:
        step_name_1 = step_1_names[0]
    else:
        step_name_1 = utils.step_name_str(yaml_stem, 0, steps_keys[0])
        step_name_1 = '___'.join(namespaces + [step_name_1])
    # NOTE: Since the names of subgraphs '*.yml' contain a period, we need to
    # escape them by enclosing the whole name in double quotes. Otherwise:
    # "Error: *.yml.gv: syntax error in line n near '.'"
        step_name_1 = f'"{step_name_1}"'

    return step_name_1


def get_workflow_outputs(args: argparse.Namespace,
                         namespaces: List[str],
                         is_root: bool,
                         yaml_stem: str,
                         steps: List[Yaml],
                         vars_workflow_output_internal: InternalOutputs,
                         graph: Graph,
                         tools: Tools,
                         step_node_name: str) -> Dict[str, Dict[str, str]]:
    # Add the outputs of each step to the workflow outputs
    workflow_outputs = {}
    steps_keys = utils.get_steps_keys(steps)
    for i, step_key in enumerate(steps_keys):
        stem = Path(step_key).stem
        tool_i = tools[stem][1]
        step_name_i = utils.step_name_str(yaml_stem, i, step_key)
        out_keys = steps[i][step_key]['out']
        for out_key in out_keys:
            # Exclude certain output files as per the comment above.
            if 'dhdl' in out_key or 'xtc' in out_key:
                continue
            out_var = f'{step_name_i}/{out_key}'
            # Avoid duplicating intermediate outputs in GraphViz
            out_key_no_namespace = out_key.split('___')[-1]
            if args.graph_show_outputs:
                case1 = (tool_i['class'] == 'Workflow') and (not out_key in [var.replace('/', '___') for var in vars_workflow_output_internal])
                # Avoid duplicating outputs from subgraphs in parent graphs.
                output_node_name = '___'.join(namespaces + [step_name_i, out_key])
                # TODO: check is_root here
                case1 = case1 and not is_root and (len(step_node_name.split('___')) + 1 == len(output_node_name.split('___')))
                case2 = (tool_i['class'] == 'CommandLineTool') and (not out_var in vars_workflow_output_internal)
                if case1 or case2:
                    (graph_gv, graph_nx) = graph
                    graph_gv.node(output_node_name, label=out_key_no_namespace, shape='box', style='rounded, filled', fillcolor='lightyellow')
                    if args.graph_label_edges:
                        graph_gv.edge(step_node_name, output_node_name, color=font_edge_color, label=out_key_no_namespace)  # Is labeling necessary?
                    else:
                        graph_gv.edge(step_node_name, output_node_name, color=font_edge_color)
                    graph_nx.add_node(output_node_name)
                    graph_nx.add_edge(step_node_name, output_node_name)
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
            workflow_outputs.update({out_name: {'type': 'File', 'outputSource': out_var}})
        #print('workflow_outputs', workflow_outputs)
    return workflow_outputs