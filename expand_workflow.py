
import argparse
import copy
import glob
from pathlib import Path
import subprocess as sub

import graphviz
import yaml


def expand_workflow(args, namespaces, dot, vars_dollar_input, tools, is_root, yaml_path):
    # Load the high-level yaml workflow file.
    yaml_stem = Path(yaml_path).stem
    with open(Path(yaml_path), 'r') as y:
        yaml_tree = yaml.load(y.read(), Loader=yaml.SafeLoader)

    steps = yaml_tree['steps']

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

    # Collect recursive dollar variable dependency mappings
    vars_dollar_tails = {}

    # Collect recursive subworkflow data
    node_1_names = []
    subdots = []

    for i, key in enumerate(steps_keys):
        stem = Path(key).stem
        # Recursively expand subworkflows, adding expanded cwl file contents to tools
        if key in subkeys:
            path = Path(key)
            if not (path.exists() and path.suffix == '.yml'):
                # TODO: Once we have defined a yml DSL schema,
                # check that the file contents actually satisfies the schema.
                raise Exception(f'Error! {path} does not exists or is not a .yml file.')
            subdot = graphviz.Digraph(name=f'cluster_{key}')
            subdot.attr(label=key)
            subworkflow_data = expand_workflow(args, namespaces + [f'{yaml_stem}_step_{i+1}_{key}'], subdot, vars_dollar_input, tools, False, path)
            subdots.append(subworkflow_data[-1])
            node_1_names.append(subworkflow_data[-2])
            # Add expanded cwl file contents to tools
            tools[stem] = (stem + '.cwl', subworkflow_data[0])

            # Initialize the above from recursive values.
            # Do not initialize inputs_workflow. See comment below.
            # inputs_workflow.update(subworkflow_data[1])
            inputs_namespaced = dict([(f'{yaml_stem}_step_{i+1}_{key}___{k}', val) for k, val in subworkflow_data[2].items()]) # _{key}_input___{k}
            inputs_file_workflow.update(inputs_namespaced)
            vars_dollar_input.update(subworkflow_data[3])
            vars_dollar_tails.update(subworkflow_data[4])
            vars_workflow_output_internal += subworkflow_data[5]

        tool_i = tools[stem][1]

        # Add run tag
        run_path = tools[stem][0]
        if steps[i][key]:
            if not 'run' in steps[i][key]:
                steps[i][key].update({'run': run_path})
        else:
            steps[i] = {key: {'run': run_path}}

        # Generate intermediate file names between steps.
        # The inputs for a given step need to come from either
        # 1. the input.yaml file for the overall workflow (extract into separate yml file) or
        # 2. the outputs from the previous steps (most recent first).
        # If there isn't an exact match, remove input_* and output_* from the
        # current step and previous steps, respectively, and then check again.
        # If there still isn't an exact match, explicit renaming may be required.

        args_provided = []
        if steps[i][key] and 'in' in steps[i][key]:
            args_provided = list(steps[i][key]['in'])
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
            steps[i][key]['in'] = dict([(key, key) for key in args_required])
        else:
            raise Exception(f'Unknown class', tool_i['class'])
        
        # Note: Some config tags are not required in the cwl files, but are in
        # fact required in the python source code! See check_mandatory_property
        # (Solution: refactor all required arguments out of config and list
        # them as explicit inputs in the cwl files, then modify the python
        # files accordingly.)
        #print(args_required)

        sub_args_provided = [arg for arg in args_required if arg in vars_dollar_tails]
        #print(sub_args_provided)

        step_name = f'{yaml_stem}_step_{i + 1}_{steps_keys[i]}' # steps_keys[i] = key
        label = key
        if args.graph_label_stepname:
            label = step_name
        if not tool_i['class'] == 'Workflow':
            dot.node(step_name, label=label, shape='box', style='rounded, filled', fillcolor='lightblue')
        elif not (steps_keys[i] in subkeys and args.graph_inline_subgraphs):
            dot.node(step_name, label=label, shape='box', style='rounded, filled', fillcolor='lightblue')

        # NOTE: sub_args_provided are handled within the args_required loop below
        for arg_key in args_provided:
            # Extract input value into separate yml file
            # Replace it here with a new variable name
            arg_val = steps[i][key]['in'][arg_key]
            in_name = f'{step_name}___{arg_key}'  # Use triple underscore for namespacing so we can split later # {step_name}_input___{arg_key}
            in_type = in_tool[arg_key]['type'].replace('?', '')  # Providing optional arguments makes them required
            if arg_val[0] == '$':
                arg_val = arg_val[1:]  # Remove $
                #print('arg_key, arg_val', arg_key, arg_val)
                if not vars_dollar_input.get(arg_val):  # if first time
                    # TODO: Add edam format info, etc.
                    inputs_workflow.update({in_name: {'type': in_type}})
                    inputs_file_workflow.update({in_name: (arg_val, in_type)})  # Add type info
                    steps[i][key]['in'][arg_key] = in_name
                    vars_dollar_input.update({arg_val: (namespaces + [step_name], arg_key)})
                else:
                    (nss, var) =  vars_dollar_input[arg_val]
                    # if recursive, add to inputs_workflow
                    if not namespaces == []:
                        inputs_workflow.update({in_name: {'type': in_type}})
                        steps[i][key]['in'][arg_key] = in_name
                        # Pass var_dollar call site info up through the recursion.
                        vars_dollar_tails.update({in_name: vars_dollar_input[arg_val]}) # {in_name, (namespaces + [step_name], var)} ?
                    else:
                        var = nss[0] + '/' + '___'.join(nss[1:] + [var])
                        steps[i][key]['in'][arg_key] = var
                    # NOTE: Do NOT add edges here! (i.e. within the subgraph)
                    # If you add an edge whose head/tail is outside of the subgraph,
                    # graphviz will segfault! Add the edges in the parent graph below.
            else:
                # TODO: Add edam format info, etc.
                inputs_workflow.update({in_name: {'type': in_type}})
                inputs_file_workflow.update({in_name: (arg_val, in_type)})  # Add type info
                steps[i][key]['in'][arg_key] = in_name
                if args.graph_show_inputs:
                    var_name = f'{step_name}/{arg_key}'
                    dot.node(var_name, label=arg_key, shape='box', style='rounded, filled', fillcolor='lightgreen')
                    dot.edge(var_name, step_name)
        
        for arg_key in args_required:
            #print('arg_key', arg_key)
            if arg_key in args_provided:
                continue  # We already covered this case above.
            match = False
            if arg_key in sub_args_provided: # Edges have been explicitly provided
                match = True
                (nss, var) = vars_dollar_tails[arg_key]
                var_slash = nss[0] + '/' + '___'.join(nss[1:] + [var])
                arg_keyval = {arg_key: var_slash}

                # NOTE: Add edges to subgraphs here (i.e. in the parent graph).
                # Otherwise, graphviz will segfault! See comment above.
                edge_node1 = nss[-1]
                edge_node2 = arg_key.split('___')[0]
                label = nss[0]

                # TODO: vars_workflow_output_internal?
            else: # We need to infer the edges
                # TODO: Figure out something better than replace. Use type and/or edam format info? Yes.
                arg_key_no_namespace = arg_key.split('___')[-1]
                arg_key_noinput = arg_key_no_namespace.replace('input_', '').replace('solute_', '').replace('energy_', 'edr_').replace('structure_', 'tpr_').replace('traj_', 'trr_')
                # gmx_trjconv_str assumes gro, but other tools assume tpr; For now,
                # instead of figuring out how to replace('structure_', 'gro_') for
                # gmx_trjconv_str only, just require the user to specify filename.
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
                            arg_keyval = {arg_key: arg_val}

                            # Determine which head and tail node to use for the new edge
                            # TODO: Check this indexing
                            edge_node1 = step_name_j
                            if tool_j['class'] == 'Workflow' and steps_keys[j] in subkeys: # i.e. if we performed recursion in step j
                                if args.graph_inline_subgraphs:
                                    edge_node1 = out_key.split('___')[0]
                            edge_node2 = step_name
                            if tool_i['class'] == 'Workflow' and steps_keys[i] in subkeys: # i.e. if we performed recursion in step i
                                if args.graph_inline_subgraphs:
                                    edge_node2 = arg_key.split('___')[0]

                            # We also need to keep track of the 'internal' output variables
                            if tool_j['class'] == 'Workflow':
                                vars_workflow_output_internal.append(out_key)
                            else:
                                vars_workflow_output_internal.append(arg_val)

                            label = out_key_no_namespace if tool_j['class'] == 'Workflow' else out_key
                    # Break out two levels, i.e. continue onto the next iteration of the outermost loop.
                            match = True
                            break
                    if match:
                        break

            if match:
                if args.graph_label_edges:
                    dot.edge(edge_node1, edge_node2, label=label)
                else:
                    dot.edge(edge_node1, edge_node2)

            if not match:
                in_name = f'{step_name}___{arg_key}'  # Use triple underscore for namespacing so we can split later # {step_name}_input___{arg_key}
                arg_keyval = {arg_key: in_name}
                # This just means we need to defer to the parent workflow.
                # There will actually be an error only if no parent supplies an
                # input value and thus there is no entry in inputs_file_workflow.
                if is_root and not in_name in inputs_file_workflow:
                    print('Error! No match found for input', i + 1, key, arg_key)

                # Add an input name to this subworkflow (only). Do not add to
                # inputs_file_workflow because this may be an internal input,
                # i.e. it may simply be split across subworkflow boundaries and
                # the input may be supplied in one of the parent workflows.
                # We also do not (necessarily) need to explicitly pass this list
                # to the parent workflow since we are expanding the 'in' tag here,
                # which should match in the parent workflow.
                in_arg = in_tool[arg_key]
                # TODO: check this
                in_type = in_arg if isinstance(in_arg, str) else in_arg['type']
                in_type = in_type.replace('?', '')  # Providing optional arguments makes them required
                inputs_workflow.update({in_name: {'type': in_type}})

            # Finally, add the inputs (for all cases).
            if steps[i][key]:
                if 'in' in steps[i][key]:
                    new_keys = dict(list(steps[i][key]['in'].items()) + list(arg_keyval.items()))
                    new_keyvals = [(k, v) if k != 'in' else (k, new_keys) for k, v in steps[i][key].items()]
                else:
                    new_keys = arg_keyval
                    new_keyvals = dict(list(steps[i][key].items()) + [('in', new_keys)])
                steps[i][key].update(new_keyvals)
            else:
                steps[i] = {key: {'in': arg_keyval}}
                
        
        # Add CommandLineTool outputs tags to workflow out tags.
        # Note: Add all output tags for now, but depending on config options,
        # not all output files will be generated. This may cause an error.
        out_keys = list(tool_i['outputs'])
        #print('out_keys', out_keys)
        if steps[i][key]:
            if 'out' in steps[i][key]:
                new_keys = steps[i][key]['out'] + out_keys
                new_keyvals = [(k, v) if k != 'out' else (k, new_keys) for k, v in steps[i][key].items()]
            else:
                new_keys = out_keys
                new_keyvals = list(steps[i][key].items()) + [('out', new_keys)]
            steps[i][key].update(new_keyvals)
        else:
            steps[i] = {key: {'out': out_keys}}

        #print()

    # Add the cluster subgraphs to the main graph, but we need to add them in
    # reverse order to trick the graphviz layout algorithm.
    for subdot in subdots[::-1]: # Reverse!
        dot.subgraph(subdot)
    # Align the cluster subgraphs using the same rank as the first node of each subgraph.
    # See https://stackoverflow.com/questions/6824431/placing-clusters-on-the-same-rank-in-graphviz
    if not node_1_names == []:
        nodes_same_rank = '\t{rank=same; ' + '; '.join(node_1_names) + '}\n'
        dot.body.append(nodes_same_rank)

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
    for i, key in enumerate(steps_keys):
        step_name = f'{yaml_stem}_step_{i + 1}_{steps_keys[i]}'
        out_keys = steps[i][key]['out']
        for out_key in out_keys:
            # Exclude certain output files as per the comment above.
            if 'dhdl' in out_key or 'xtc' in out_key:
                continue
            out_var = f'{step_name}/{out_key}'
            out_name = f'{step_name}___{out_key}'  # Use triple underscore for namespacing so we can split later
            # Avoid duplicating intermediate outputs in GraphViz
            out_key_no_namespace = out_key.split('___')[-1]
            if args.graph_show_outputs:
                case1 = (tool_i['class'] == 'Workflow') and (not out_key in [var.replace('/', '___') for var in vars_workflow_output_internal])
                # Avoid duplicating outputs from subgraphs in the main graph.
                case1 = case1 and not args.graph_inline_subgraphs and not is_root
                case2 = (tool_i['class'] == 'CommandLineTool') and (not out_var in vars_workflow_output_internal)
                if case1 or case2:
                    dot.node(out_name, label=out_key, shape='box', style='rounded, filled', fillcolor='lightyellow')
                    if args.graph_label_edges:
                        dot.edge(step_name, out_name, label=out_key_no_namespace)  # Is labeling necessary?
                    else:
                        dot.edge(step_name, out_name)
            # Exclude intermediate 'output' files.
            if out_var in vars_workflow_output_internal and not args.cwl_output_intermediate_files:
                continue
            #print('out_name', out_name)
            outputs_workflow.update({out_name: {'type': 'File', 'outputSource': out_var}})
        #print('outputs_workflow', outputs_workflow)
    yaml_tree.update({'outputs': outputs_workflow})

    # Finally, rename the steps to be unique by prepending step_*_
    # and convert the list of steps into a dict
    steps_dict = {}
    for i, key in enumerate(steps_keys):
        step_name = f'{yaml_stem}_step_{i + 1}_{steps_keys[i]}'
        #steps[i] = {step_name: steps[i][key]}
        steps_dict.update({step_name: steps[i][key]})
    yaml_tree.update({'steps': steps_dict})
    
    # Dump the workflow inputs to a separate yml file.
    temp = {}
    for key, (val, in_type) in inputs_file_workflow.items():
        new_keyval = {key: val}  # if string
        if 'File' in in_type:
            new_keyval = {key: {'class': 'File', 'path': val, 'format': 'https://edamontology.org/format_2330'}}
        temp.update(new_keyval)
            
    dump_options = {'line_break': '\n', 'indent': 2}
    yaml_content = yaml.dump(temp, sort_keys=False, **dump_options)
    with open(f'{yaml_stem}_inputs.yml', 'w') as inp:
        inp.write(yaml_content)

    # Dump the expanded yaml file to disk.
    # Use sort_keys=False to preserve the order of the steps.
    dump_options = {'line_break': '\n', 'indent': 2}
    yaml_content = yaml.dump(yaml_tree, sort_keys=False, **dump_options)
    with open(f'{yaml_stem}.cwl', 'w') as w:
        w.write('#!/usr/bin/env cwl-runner\n')
        w.write(''.join(yaml_content))
    
    print(f'Finished expanding {yaml_path}')
    
    if args.cwl_validate:
        print(f'Validating {yaml_stem}.cwl ...')
        cmd = ['cwltool', '--validate', f'{yaml_stem}.cwl']
        sub.run(cmd)
    
    # Note: We do not necessarily need to return inputs_workflow.
    # 'Internal' inputs are encoded in yaml_tree. See Comment above.
    node_1_name = f'{yaml_stem}_step_1_{steps_keys[0]}'  # Solely for graphviz layout purposes.
    return (yaml_tree, inputs_workflow, inputs_file_workflow, vars_dollar_input, vars_dollar_tails, vars_workflow_output_internal, is_root, node_1_name, dot)


def main():
    parser = argparse.ArgumentParser(prog='main', description='Convert a high-level yaml workflow file to CWL.')
    parser.add_argument('--yaml', type=str, required=True,
                        help='Yaml workflow file')
    parser.add_argument('--cwl_dir', type=str, required=False, default='biobb',
                        help='Directory which contains the CWL CommandLineTools and/or Workflows')
    parser.add_argument('--cwl_output_intermediate_files', type=bool, required=False, default=False,
                        help='Enable output files which are used between steps (for debugging).')
    parser.add_argument('--cwl_run', type=bool, required=False, default=False,
                        help='After generating the cwl file, run it.')
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
    parser.add_argument('--graph_inline_subgraphs', type=bool, required=False, default=False,
                        help='Controls whether subgraphs are displayed separately or positioned within the main graph.')
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
    for cwl_path in cwl_paths_sorted:
        #print(cwl_path)
        try:
            with open(cwl_path, 'r') as f:
              tool = yaml.load(f.read(), Loader=yaml.SafeLoader)
            stem = Path(cwl_path).stem
            #print(stem)
            tools_cwl[stem] = (cwl_path, tool)
            #print(tool)
        except yaml.scanner.ScannerError as se:
            pass
            # There are two cwl files that throw this error, but they are both legacy, so...
            #print(cwl_path)
            #print(se)

    # Collect the explicit $ internal workflow input variables
    vars_dollar_input = {}

    rootdot = graphviz.Digraph(name=args.yaml)
    rootdot.attr(newrank='True') # See graphviz layout comment above.
    with rootdot.subgraph(name=f'cluster_{args.yaml}') as subdot:
        subdot.attr(label=args.yaml)
        workflow_data = expand_workflow(args, [], subdot, vars_dollar_input, tools_cwl, True, args.yaml)
    # Render the GraphViz diagram
    rootdot.render()
    #rootdot.view() # viewing does not work on headless machines (and requires xdg-utils)
    
    if args.cwl_run:
        yaml_stem = Path(args.yaml).stem
        print(f'Running {yaml_stem}.cwl ...')
        cmd = ['cwltool', '--cachedir', 'cachedir','--outdir', 'outdir', f'{yaml_stem}.cwl', f'{yaml_stem}_inputs.yml']
        sub.run(cmd)


if __name__ == '__main__':
    main()
