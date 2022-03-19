
import argparse
import glob
from pathlib import Path
import requests
import subprocess as sub

import graphviz
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


def perform_edge_inference(args, tools, steps_keys, subkeys, yaml_stem, i, steps_i, arg_key, graph, is_root, namespaces, inputs_workflow, inputs_file_workflow, vars_workflow_output_internal):
    match = False
    # TODO: Figure out something better than replace. Use type and/or edam format info? Yes.
    arg_key_no_namespace = arg_key.split('___')[-1]
    arg_key_noinput = arg_key_no_namespace.replace('input_', '').replace('solute_', '').replace('energy_', 'edr_').replace('structure_', 'tpr_').replace('traj_', 'trr_')
    # gmx_trjconv_str assumes gro, but other tools assume tpr; For now,
    # instead of figuring out how to replace('structure_', 'gro_') for
    # gmx_trjconv_str only, just require the user to specify filename.
    step_key = steps_keys[i]
    tool_i = tools[Path(step_key).stem][1]
    step_name_i = f'{yaml_stem}_step_{i + 1}_{step_key}'
    for j in range(0, i)[::-1]:  # Reverse order!
        tool_j = tools[Path(steps_keys[j]).stem][1]
        out_keys = list(tool_j['outputs'])[::-1] # Reverse order!
        for out_key in out_keys:
            out_key_no_namespace = out_key.split('___')[-1]
            if arg_key_noinput == out_key_no_namespace.replace('output_', ''):
                #print('match!', j)  # We found a match!
                # Generate a new namespace for out_key using the step number and add to inputs
                step_name_j = f'{yaml_stem}_step_{j + 1}_{steps_keys[j]}'
                arg_val = f'{step_name_j}/{out_key}'

                # We also need to keep track of the 'internal' output variables
                if tool_j['class'] == 'Workflow':
                    vars_workflow_output_internal.append(out_key)
                else:
                    vars_workflow_output_internal.append(arg_val)

                # Determine which head and tail node to use for the new edge
                # First we need to extract the embedded namespaces
                nss_embedded1 = out_key.split('___')[:-1]
                nss_embedded2 = arg_key.split('___')[:-1]
                nss1 = namespaces + [step_name_j] + nss_embedded1
                nss2 = namespaces + [step_name_i] + nss_embedded2
                # TODO: check this
                label = out_key_no_namespace if tool_j['class'] == 'Workflow' else out_key
                add_graph_edge(args, graph, nss1, nss2, label)

                arg_keyval = {arg_key: arg_val}
                steps_i = add_yamldict_keyval(steps_i, step_key, 'in', arg_keyval)
                return steps_i  # Short circuit

    if not match:
        in_name = f'{step_name_i}___{arg_key}'  # Use triple underscore for namespacing so we can split later # {step_name_i}_input___{arg_key}

        # This just means we need to defer to the parent workflow.
        # There will actually be an error only if no parent supplies an
        # input value and thus there is no entry in inputs_file_workflow.
        if is_root and not in_name in inputs_file_workflow:
            print('Error! No match found for input', i + 1, step_key, arg_key)

        # Add an input name to this subworkflow (only). Do not add to
        # inputs_file_workflow because this may be an internal input,
        # i.e. it may simply be split across subworkflow boundaries and
        # the input may be supplied in one of the parent workflows.
        # We also do not (necessarily) need to explicitly pass this list
        # to the parent workflow since we are expanding the 'in' tag here,
        # which should match in the parent workflow.
        in_tool = tool_i['inputs']
        in_arg = in_tool[arg_key]
        # TODO: check this
        in_type = in_arg if isinstance(in_arg, str) else in_arg['type']
        in_type = in_type.replace('?', '')  # Providing optional arguments makes them required
        inputs_workflow.update({in_name: {'type': in_type}})

        arg_keyval = {arg_key: in_name}
        steps_i = add_yamldict_keyval(steps_i, step_key, 'in', arg_keyval)
        return steps_i


def expand_workflow(args, namespaces, subgraphs, vars_dollar_defs, tools, is_root, yaml_path, yml_paths):
    # Load the high-level yaml workflow file.
    yaml_stem = Path(yaml_path).stem
    with open(Path(yaml_path), 'r') as y:
        yaml_tree = yaml.safe_load(y.read())

    steps = get_steps(yaml_tree, yaml_path)

    # Get the dictionary key (i.e. the name) of each step.
    steps_keys = []
    for step in steps:
        steps_keys += list(step)
    #print(steps_keys)

    subkeys = [key for key in steps_keys if key not in tools]

    # Add headers
    yaml_tree['cwlVersion'] = 'v1.0'
    yaml_tree['class'] = 'Workflow'

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

    # Collect workflow input parameters
    inputs_workflow = {}
    inputs_file_workflow = {}
    
    # Collect the internal workflow output variables
    vars_workflow_output_internal = []

    # Collect recursive dollar variable call sites.
    vars_dollar_calls = {}

    # Collect recursive subworkflow data
    step_1_names = []
    sibling_subgraphs = []

    # Collect labshare plugin_ids
    plugin_ids = []

    graph = subgraphs[-1] # Get the current graph

    for i, step_key in enumerate(steps_keys):
        stem = Path(step_key).stem
        # Recursively expand subworkflows, adding expanded cwl file contents to tools
        if step_key in subkeys:
            #path = Path(step_key)
            path = Path(yml_paths[Path(step_key).stem])
            if not (path.exists() and path.suffix == '.yml'):
                # TODO: Once we have defined a yml DSL schema,
                # check that the file contents actually satisfies the schema.
                raise Exception(f'Error! {path} does not exist or is not a .yml file.')
            subgraph = graphviz.Digraph(name=f'cluster_{step_key}')
            subgraph.attr(label=step_key) # str(path)
            subgraph.attr(color='lightblue')  # color of outline
            subworkflow_data = expand_workflow(args, namespaces + [f'{yaml_stem}_step_{i+1}_{step_key}'], subgraphs + [subgraph], vars_dollar_defs, tools, False, path, yml_paths)
            sibling_subgraphs.append(subworkflow_data[-1]) # TODO: Just subgraph?
            step_1_names.append(subworkflow_data[-2])
            # Add expanded cwl file contents to tools
            tools[stem] = (stem + '.cwl', subworkflow_data[0])

            # Initialize the above from recursive values.
            # Do not initialize inputs_workflow. See comment below.
            # inputs_workflow.update(subworkflow_data[1])
            inputs_namespaced = dict([(f'{yaml_stem}_step_{i+1}_{step_key}___{k}', val) for k, val in subworkflow_data[2].items()]) # _{step_key}_input___{k}
            inputs_file_workflow.update(inputs_namespaced)
            vars_dollar_defs.update(subworkflow_data[3])
            vars_dollar_calls.update(subworkflow_data[4])
            vars_workflow_output_internal += subworkflow_data[5]

        tool_i = tools[stem][1]

        plugin_id = None
        if args.cwl_run_slurm:
            if not (step_key in subkeys):
                # i.e. If this is either a primitive CommandLineTool and/or
                # a 'primitive' Workflow that we did NOT recursively generate.
                plugin_id = upload_labshare_plugin(args.compute_url, tool_i, stem)
            else:
                plugin_id = subworkflow_data[6]
            plugin_ids.append(plugin_id)

        # Add run tag
        run_path = tools[stem][0]
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

        step_name_i = f'{yaml_stem}_step_{i + 1}_{step_key}'
        label = step_key
        if args.graph_label_stepname:
            label = step_name_i
        step_node_name = '___'.join(namespaces + [step_name_i])
        if not tool_i['class'] == 'Workflow':
            graph.node(step_node_name, label=label, shape='box', style='rounded, filled', fillcolor='lightblue')
        elif not (step_key in subkeys and len(namespaces) < args.graph_inline_depth):
            nssnode = namespaces + [step_name_i]
            nssnode = nssnode[:(1 + args.graph_inline_depth)]
            step_node_name = '___'.join(nssnode)
            graph.node(step_node_name, label=label, shape='box', style='rounded, filled', fillcolor='lightblue')

        # NOTE: sub_args_provided are handled within the args_required loop below
        for arg_key in args_provided:
            # Extract input value into separate yml file
            # Replace it here with a new variable name
            arg_val = steps[i][step_key]['in'][arg_key]
            in_name = f'{step_name_i}___{arg_key}'  # Use triple underscore for namespacing so we can split later # {step_name_i}_input___{arg_key}
            in_type = in_tool[arg_key]['type'].replace('?', '')  # Providing optional arguments makes them required
            if arg_val[0] == '$':
                arg_val = arg_val[1:]  # Remove $
                #print('arg_key, arg_val', arg_key, arg_val)
                # NOTE: There can only be one definition, but multiple call sites.
                if not vars_dollar_defs.get(arg_val):
                    # if first time encountering arg_val, i.e. if defining
                    # TODO: Add edam format info, etc.
                    inputs_workflow.update({in_name: {'type': in_type}})
                    inputs_file_workflow.update({in_name: (arg_val, in_type)})  # Add type info
                    steps[i][step_key]['in'][arg_key] = in_name
                    vars_dollar_defs.update({arg_val: (namespaces + [step_name_i], arg_key)})
                    # TODO: Show input node?
                else:
                    (nss_def_init, var) =  vars_dollar_defs[arg_val]

                    nss_def_embedded = var.split('___')[:-1]
                    nss_call_embedded = arg_key.split('___')[:-1]
                    nss_def = nss_def_init + nss_def_embedded
                    # [step_name_i] is correct; nss_def_init already contains step_name_j from the recursive call
                    nss_call = namespaces + [step_name_i] + nss_call_embedded

                    nss_def_inits, nss_def_tails = partition_by_lowest_common_ancestor(nss_def, nss_call)
                    nss_call_inits, nss_call_tails = partition_by_lowest_common_ancestor(nss_call, nss_def)
                    assert nss_def_inits == nss_call_inits
                    
                    # If the call site of an explicit edge is at the same level
                    # or deeper into the recursion than the definition, or if
                    # they are in completely different branches, then we
                    # need to create inputs in all of the intervening cwl files
                    # so we can pass in the values from the outer scope(s). Here,
                    # we simply need to use in_name and add to inputs_workflow
                    # and vars_dollar_calls. The outer scope(s) are handled by
                    # the `if not match` clause in perform_edge_inference().
                    if len(nss_call_tails) >= len(nss_def_tails) or nss_call_inits == []:
                        inputs_workflow.update({in_name: {'type': in_type}})
                        steps[i][step_key]['in'][arg_key] = in_name
                        # Store var_dollar call site info up through the recursion.
                        vars_dollar_calls.update({in_name: vars_dollar_defs[arg_val]}) # {in_name, (namespaces + [step_name_i], var)} ?
                    else:
                        var_slash = nss_def_tails[0] + '/' + '___'.join(nss_def_tails[1:] + [var])
                        steps[i][step_key]['in'][arg_key] = var_slash

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
                    add_graph_edge(args, graph_init, nss_def, nss_call, label)
            else:
                # TODO: Add edam format info, etc.
                inputs_workflow.update({in_name: {'type': in_type}})
                inputs_file_workflow.update({in_name: (arg_val, in_type)})  # Add type info
                steps[i][step_key]['in'][arg_key] = in_name
                if args.graph_show_inputs:
                    input_node_name = '___'.join(namespaces + [step_name_i, arg_key])
                    graph.node(input_node_name, label=arg_key, shape='box', style='rounded, filled', fillcolor='lightgreen')
                    graph.edge(input_node_name, step_node_name, color=font_edge_color)
        
        for arg_key in args_required:
            #print('arg_key', arg_key)
            if arg_key in args_provided:
                continue  # We already covered this case above.
            if arg_key in sub_args_provided: # Edges have been explicitly provided
                # We have now returned from the recursion and need to pass values in.
                # Extract the stored defs namespaces from vars_dollar_calls.
                # (See massive comment above.)
                (nss_def_init, var) = vars_dollar_calls[arg_key]

                nss_def_embedded = var.split('___')[:-1]
                nss_call_embedded = arg_key.split('___')[:-1]
                nss_def = nss_def_init + nss_def_embedded
                # [step_name_i] is correct; nss_def_init already contains step_name_j from the recursive call
                nss_call = namespaces + [step_name_i] + nss_call_embedded

                nss_def_inits, nss_def_tails = partition_by_lowest_common_ancestor(nss_def, nss_call)
                nss_call_inits, nss_call_tails = partition_by_lowest_common_ancestor(nss_call, nss_def)
                assert nss_def_inits == nss_call_inits

                var_slash = nss_def_tails[0] + '/' + '___'.join(nss_def_tails[1:] + [var])
                arg_keyval = {arg_key: var_slash}
                steps[i] = add_yamldict_keyval(steps[i], step_key, 'in', arg_keyval)

                # NOTE: We already added an edge to the appropriate subgraph above.
                # TODO: vars_workflow_output_internal?
            else:
                steps[i] = perform_edge_inference(args, tools, steps_keys, subkeys, yaml_stem, i, steps[i], arg_key, graph, is_root, namespaces, inputs_workflow, inputs_file_workflow, vars_workflow_output_internal)
        
        # Add CommandLineTool outputs tags to workflow out tags.
        # Note: Add all output tags for now, but depending on config options,
        # not all output files will be generated. This may cause an error.
        out_keys = list(tool_i['outputs'])
        #print('out_keys', out_keys)
        steps[i] = add_yamldict_keyval(steps[i], step_key, 'out', out_keys)

        #print()

    # Add the cluster subgraphs to the main graph, but we need to add them in
    # reverse order to trick the graphviz layout algorithm.
    if len(namespaces) < args.graph_inline_depth:
        for sibling in sibling_subgraphs[::-1]: # Reverse!
            graph.subgraph(sibling)
    # Align the cluster subgraphs using the same rank as the first node of each subgraph.
    # See https://stackoverflow.com/questions/6824431/placing-clusters-on-the-same-rank-in-graphviz
    if len(namespaces) < args.graph_inline_depth:
        step_1_names_display = [name for name in step_1_names if len(name.split('___')) < 2 + args.graph_inline_depth]
        if len(step_1_names_display) > 1:
            nodes_same_rank = '\t{rank=same; ' + '; '.join(step_1_names_display) + '}\n'
            graph.body.append(nodes_same_rank)
    if steps_keys[0] in subkeys:
        step_name_1 = step_1_names[0]
    else:
        step_name_1 = f'{yaml_stem}_step_1_{steps_keys[0]}'
        step_name_1 = '___'.join(namespaces + [step_name_1])
    # NOTE: Since the names of subgraphs '*.yml' contain a period, we need to
    # escape them by enclosing the whole name in double quotes. Otherwise:
    # "Error: *.yml.gv: syntax error in line n near '.'"
        step_name_1 = f'"{step_name_1}"'

    # Add the provided inputs of each step to the workflow inputs
    temp = {}
    for k, v in inputs_workflow.items():
        # TODO: Remove this heuristic after we add the format info above.
        new_keyval = {k: v['type']}
        if 'mdin' in k and 'File' in v['type']:
            newval = dict(list(v.items()) + list({'format': 'https://edamontology.org/format_2330'}.items()))
            new_keyval = {k: newval}
        temp.update(new_keyval)
    yaml_tree.update({'inputs': temp})

    # Add the outputs of each step to the workflow outputs
    vars_workflow_output_internal = list(set(vars_workflow_output_internal))  # Get uniques
    outputs_workflow = {}
    for i, step_key in enumerate(steps_keys):
        step_name_i = f'{yaml_stem}_step_{i + 1}_{steps_keys[i]}'
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
                case1 = case1 and not is_root and (len(step_node_name.split('___')) + 1 == len(output_node_name.split('___')))
                case2 = (tool_i['class'] == 'CommandLineTool') and (not out_var in vars_workflow_output_internal)
                if case1 or case2:
                    graph.node(output_node_name, label=out_key_no_namespace, shape='box', style='rounded, filled', fillcolor='lightyellow')
                    if args.graph_label_edges:
                        graph.edge(step_node_name, output_node_name, color=font_edge_color, label=out_key_no_namespace)  # Is labeling necessary?
                    else:
                        graph.edge(step_node_name, output_node_name, color=font_edge_color)
            # Exclude intermediate 'output' files.
            if out_var in vars_workflow_output_internal and not args.cwl_output_intermediate_files:
                continue
            out_name = f'{step_name_i}___{out_key}'  # Use triple underscore for namespacing so we can split later
            #print('out_name', out_name)
            outputs_workflow.update({out_name: {'type': 'File', 'outputSource': out_var}})
        #print('outputs_workflow', outputs_workflow)
    yaml_tree.update({'outputs': outputs_workflow})

    # Finally, rename the steps to be unique by prepending step_*_
    # and convert the list of steps into a dict
    steps_dict = {}
    for i, step_key in enumerate(steps_keys):
        step_name_i = f'{yaml_stem}_step_{i + 1}_{step_key}'
        #steps[i] = {step_name_i: steps[i][step_key]}
        steps_dict.update({step_name_i: steps[i][step_key]})
    yaml_tree.update({'steps': steps_dict})

    # Dump the workflow inputs to a separate yml file.
    yaml_inputs = {}
    for key, (val, in_type) in inputs_file_workflow.items():
        new_keyval = {key: val}  # if string
        if 'File' in in_type:
            new_keyval = {key: {'class': 'File', 'path': val, 'format': 'https://edamontology.org/format_2330'}}
        yaml_inputs.update(new_keyval)
            
    dump_options = {'line_break': '\n', 'indent': 2}
    yaml_content = yaml.dump(yaml_inputs, sort_keys=False, **dump_options)
    with open(f'{yaml_stem}_inputs.yml', 'w') as inp:
        inp.write(yaml_content)

    # Dump the expanded yaml file to disk.
    # Use sort_keys=False to preserve the order of the steps.
    dump_options = {'line_break': '\n', 'indent': 2}
    yaml_content = yaml.dump(yaml_tree, sort_keys=False, **dump_options)
    with open(f'{yaml_stem}.cwl', 'w') as w:
        w.write('#!/usr/bin/env cwl-runner\n')
        w.write(''.join(yaml_content))

    if args.cwl_run_slurm:
        # Convert the expanded yaml file to json for labshare Compute.
        # Replace 'run' with plugin:id
        import copy
        yaml_tree_run = copy.deepcopy(yaml_tree)
        for i, step_key in enumerate(steps_keys):
            step_name_i = f'{yaml_stem}_step_{i + 1}_{step_key}'
            run_val = f'plugin:{plugin_ids[i]}'
            yaml_tree_run['steps'][step_name_i]['run'] = run_val
        if key in subkeys: # and not is_root, but the former implies the latter
            plugin_id = upload_labshare_plugin(args.compute_url, yaml_tree_run, stem)
        if is_root:
            compute_workflow = {
                "driver": "slurm",
                "name": yaml_stem,
                "cwlJobInputs": yaml_inputs,
                **yaml_tree_run
            }
            # Use http POST request to upload a primitive CommandLineTool / define a plugin and get its id hash.
            response = requests.post(args.compute_url + '/compute/workflows', json = compute_workflow)
            r_json = response.json()
            print('post response')
            print(r_json)
            if 'id' not in r_json:
                raise Exception(f'Error! Labshare workflow upload failed for {yaml_stem}.')
            print_labshare_plugins(args.compute_url)

    print(f'Finished expanding {yaml_path}')
    
    if args.cwl_validate:
        print(f'Validating {yaml_stem}.cwl ...')
        cmd = ['cwltool', '--validate', f'{yaml_stem}.cwl']
        sub.run(cmd)
    
    # Note: We do not necessarily need to return inputs_workflow.
    # 'Internal' inputs are encoded in yaml_tree. See Comment above.
    return (yaml_tree, inputs_workflow, inputs_file_workflow, vars_dollar_defs, vars_dollar_calls, vars_workflow_output_internal, plugin_id, is_root, step_name_1, graph)


def main():
    parser = argparse.ArgumentParser(prog='main', description='Convert a high-level yaml workflow file to CWL.')
    parser.add_argument('--yaml', type=str, required=True,
                        help='Yaml workflow file')
    parser.add_argument('--cwl_dir', type=str, required=False, default='biobb',
                        help='Directory which contains the CWL CommandLineTools and/or Workflows')
    parser.add_argument('--cwl_output_intermediate_files', type=bool, required=False, default=False,
                        help='Enable output files which are used between steps (for debugging).')
    parser.add_argument('--cwl_run_local', type=bool, required=False, default=False,
                        help='After generating the cwl file(s), run it on localhost.')
    # NOTE: If cwl_run_slurm is enabled, you MUST enable cwl_inline_subworkflows!
    # Plugins with 'class: Workflow' (i.e. subworkflows) are not currently supported.
    parser.add_argument('--cwl_run_slurm', type=bool, required=False, default=False,
                        help='After generating the cwl file, run it on labshare using the slurm driver.')
    parser.add_argument('--compute_url', type=bool, required=False, default=False,
                        help='The URL associated with the labshare slurm driver.')
    parser.add_argument('--cwl_inline_subworkflows', type=bool, required=False, default=False,
                        help='Before generating the cwl file, inline all subworkflows.')
    parser.add_argument('--cwl_validate', type=bool, required=False, default=False,
                        help='After generating the cwl file, validate it.')
    parser.add_argument('--graph_label_edges', type=bool, required=False, default=False,
                        help='Label the graph edges with the name of the intermediate input/output.')
    parser.add_argument('--graph_label_stepname', type=bool, required=False, default=False,
                        help='Prepend the step name to each step node.')
    parser.add_argument('--graph_show_inputs', type=bool, required=False, default=False,
                        help='Add nodes to the graph representing the workflow inputs.')
    parser.add_argument('--graph_show_outputs', type=bool, required=False, default=False,
                        help='Add nodes to the graph representing the workflow outputs.')
    parser.add_argument('--graph_inline_depth', type=int, required=False, default=100,  # 100 == all
                        help='Controls the depth of subgraphs which are displayed separately or positioned within the main graph.')
    args = parser.parse_args()

    # Load ALL of the tools.
    tools_cwl = {}
    pattern_cwl = str(Path(args.cwl_dir) / '**/*.cwl')
    #print(pattern_cwl)
    # Note that there is a current and a legacy copy of each cwl file for each tool.
    # The only difference appears to be that some legacy parameters are named 
    # *_file as opposed to *_path. Since glob does NOT return the results in
    # any particular order, and since we are using stem as our dict key, current
    # files may be overwritten with legacy files (and vice versa), resulting in
    # an inconsistent naming scheme. Since legacy files are stored in an additional
    # subdirctory, if we sort the paths by descending length, we can overwrite
    # the dict entries of the legacy files.
    cwl_paths_sorted = sorted(glob.glob(pattern_cwl, recursive=True), key=len, reverse=True)

    # Delete plugins previously uploaded to labshare.
    if args.cwl_run_slurm and Path('plugin_ids').exists():
        with open('plugin_ids', 'r') as f:
            ids = f.read().splitlines()
        for id in ids:
            response = requests.delete(args.compute_url + '/compute/plugins/' + id)
        sub.run(['rm', 'plugin_ids'])

    for cwl_path in cwl_paths_sorted:
        #print(cwl_path)
        try:
            with open(cwl_path, 'r') as f:
              tool = yaml.safe_load(f.read())
            stem = Path(cwl_path).stem
            #print(stem)
            # Add / overwrite stdout and stderr
            tool.update({'stdout': f'{stem}.out'})
            tool.update({'stderr': f'{stem}.err'})
            tools_cwl[stem] = (cwl_path, tool)
            #print(tool)
        except yaml.scanner.ScannerError as se:
            pass
            # There are two cwl files that throw this error, but they are both legacy, so...
            #print(cwl_path)
            #print(se)

    # Glob all of the yml files too, so we don't have to deal with relative paths.
    pattern_yml = str(Path(args.cwl_dir) / '**/*.yml')
    yml_paths_sorted = sorted(glob.glob(pattern_yml, recursive=True), key=len, reverse=True)
    yml_paths = {}
    for yml_path in yml_paths_sorted:
        stem = Path(yml_path).stem
        yml_paths[stem] = yml_path

    # Collect the explicit $ internal workflow input variables
    vars_dollar_defs = {}

    yaml_path = args.yaml
    if args.cwl_inline_subworkflows:
        steps_inlined = inline_sub_steps(yaml_path, tools_cwl, yml_paths)
        with open(Path(yaml_path), 'r') as y:
            yaml_tree = yaml.safe_load(y.read())
        yaml_tree['steps'] = steps_inlined
        dump_options = {'line_break': '\n', 'indent': 2}
        yaml_content = yaml.dump(yaml_tree, sort_keys=False, **dump_options)
        yaml_path = Path(args.yaml).stem + '_inline.yml'
        with open(yaml_path, 'w') as y:
            y.write(yaml_content)

    rootgraph = graphviz.Digraph(name=yaml_path)
    rootgraph.attr(newrank='True') # See graphviz layout comment above.
    rootgraph.attr(bgcolor="transparent") # Useful for making slides
    rootgraph.attr(fontcolor=font_edge_color)
    #rootgraph.attr(rankdir='LR') # When --graph_inline_depth 1, this usually looks better.
    with rootgraph.subgraph(name=f'cluster_{yaml_path}') as subgraph:
        subgraph.attr(label=yaml_path)
        subgraph.attr(color='lightblue')  # color of cluster subgraph outline
        workflow_data = expand_workflow(args, [], [subgraph], vars_dollar_defs, tools_cwl, True, yaml_path, yml_paths)
    # Render the GraphViz diagram
    rootgraph.render(format='png') # Default pdf. See https://graphviz.org/docs/outputs/
    #rootgraph.view() # viewing does not work on headless machines (and requires xdg-utils)
    
    if args.cwl_run_local:
        yaml_stem = Path(args.yaml).stem
        yaml_stem = yaml_stem + '_inline' if args.cwl_inline_subworkflows else yaml_stem
        print(f'Running {yaml_stem}.cwl ...')
        cmd = ['cwltool', '--cachedir', 'cachedir','--outdir', 'outdir', f'{yaml_stem}.cwl', f'{yaml_stem}_inputs.yml']
        sub.run(cmd)


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


def get_steps(yaml_tree, yaml_path):
    backend = None
    if 'backends' in yaml_tree:
        if 'default_backend' in yaml_tree:
            backend = yaml_tree['default_backend']
        if backend is None:
            raise Exception(f'Error! No backend in {yaml_path}!')
        if backend not in yaml_tree['backends']:
            raise Exception(f'Error! No steps for backend {backend} in {yaml_path}!')
        steps = yaml_tree['backends'][backend]['steps']
    elif 'steps' in yaml_tree:
        steps = yaml_tree['steps']
    else:
        raise Exception(f'Error! No backends and/or steps in {yaml_path}!')
    return steps


def inline_sub_steps(yaml_path, tools, yml_paths):
    # Load the high-level yaml workflow file.
    with open(Path(yaml_path), 'r') as y:
        yaml_tree = yaml.safe_load(y.read())

    steps = get_steps(yaml_tree, yaml_path)

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


def upload_labshare_plugin(compute_url, tool, stem):
        # Convert the expanded yaml file to json for labshare Compute.
        # First remove $ in $namespaces and $schemas (anywhere else?)
        dump_options = {'line_break': '\n', 'indent': 2}
        tool_no_dollar = yaml.dump(tool, sort_keys=False, **dump_options).replace('$namespaces', 'namespaces').replace('$schemas', 'schemas')
        compute_plugin = {
            # Add unique 'id' below
            'cwlScript': yaml.safe_load(tool_no_dollar)  # This effectively copies tool
        }

        # Use http POST request to upload a primitive CommandLineTool / define a plugin and get its id hash.
        response = requests.post(compute_url + '/compute/plugins', json = compute_plugin)
        r_json = response.json()
        if 'id' not in r_json:
            print('post response')
            print(r_json)
            raise Exception(f'Error! Labshare plugin upload failed for {stem}.')

        compute_plugin.update({'id': r_json['id']}) # Necessary ?
        # Save the plugin ids so we can delete them the next time we enable args.cwl_run_slurm
        with open('plugin_ids', 'a') as f:
            f.write(r_json['id'] + '\n')
        return r_json['id']


def print_labshare_plugins(compute_url):
    r = requests.get(compute_url + '/compute/plugins/')
    for j in r.json():
        print(j)
    print(len(r.json()))


if __name__ == '__main__':
    main()