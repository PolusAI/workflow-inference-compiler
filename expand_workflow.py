
import argparse
import copy
import glob
from pathlib import Path
import subprocess as sub

import graphviz
import yaml


def expand_workflow(args, dot, tools_cwl, is_root, yaml_path):
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
    
    # Recursively expand subworkflows, if any.
    subkeys = [key for key in steps_keys if key not in tools_cwl]
    tools_yml = {}
    subworkflows = {}
    for key in subkeys:
        path = Path(key)
        stem_no_ext = str(path.with_suffix(''))
        if not (path.exists() and path.suffix == '.yml'):
            # TODO: Once we have defined a yml DSL schema,
            # check that the file contents actually satisfies the schema.
            raise Exception(f'Error! {path} does not exists or is not a .yml file.')
        with dot.subgraph(name=f'cluster_{key}') as subdot:
            subdot.attr(label=stem_no_ext)
            subworkflow_data = expand_workflow(args, subdot, tools_cwl, False, path)
        (sub_yaml_tree, sub_inputs_file_workflow, sub_vars_workflow_output_internal, sub_is_root, sub_dot) = subworkflow_data
        stem = stem_no_ext + '.cwl'
        tools_yml[key] = (stem, sub_yaml_tree)
        subworkflows[key] = subworkflow_data

    tools_cwl_yml = dict(list(tools_cwl.items()) + list(tools_yml.items()))
    # Restrict to the subset of tools that we are actually using.
    run_paths = [tools_cwl_yml[key][0] for key in steps_keys]
    tools = [tools_cwl_yml[key][1] for key in steps_keys]

    # Add headers
    yaml_tree['cwlVersion'] = 'v1.0'
    yaml_tree['class'] = 'Workflow'
    # If there is at least one subworkflow, add a SubworkflowFeatureRequirement
    if any([tool['class'] == 'Workflow' for tool in tools]):
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

    # Collect the internal workflow input variables
    vars_workflow_input = {}
    
    # Collect the internal workflow output variables
    vars_workflow_output_internal = []

    for i, key in enumerate(steps_keys):
        # Initialize the above from recursive values.
        if key in subkeys:
            sub_namespaced = dict([(f'step_{i+1}_{key}_input___{k}', val) for k, val in subworkflows[key][1].items()])
            inputs_file_workflow.update(sub_namespaced)
            vars_workflow_output_internal += subworkflows[key][2]
        # Add run tag
        if steps[i][key]:
            if not 'run' in steps[i][key]:
                steps[i][key].update({'run': run_paths[i]})
        else:
            steps[i] = {key: {'run': run_paths[i]}}

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

        in_tool = tools[i]['inputs']
        #print(list(in_tool.keys()))
        if tools[i]['class'] == 'CommandLineTool':
            args_required = [arg for arg in in_tool if not (in_tool[arg].get('default') or in_tool[arg]['type'][-1] == '?')]
        elif tools[i]['class'] == 'Workflow':
            args_required = [arg for arg in in_tool]
            
            # Add the inputs. For now, assume that all sub-workflows have been
            # auto-generated from a previous application of this script, and
            # thus that all of their transitive inputs have been satisfied.
            # (i.e. simply combine the input yml files using the cat command.)
            steps[i][key]['in'] = dict([(key, key) for key in args_required])
        else:
            raise Exception(f'Unknown class', tools[i]['class'])
        
        # Note: Some config tags are not required in the cwl files, but are in
        # fact required in the python source code! See check_mandatory_property
        # (Solution: refactor all required arguments out of config and list
        # them as explicit inputs in the cwl files, then modify the python
        # files accordingly.)
        #print(args_required)
        
        step_name = f'step_{i + 1}_{steps_keys[i]}'
        if not tools[i]['class'] == 'Workflow':
            dot.node(step_name, shape='box', style='rounded, filled', fillcolor='lightblue')
        elif not steps_keys[i] in subkeys or (not args.graph_inline_subgraphs and args.graph_show_outputs):
            dot.node(step_name, shape='box', style='rounded, filled', fillcolor='lightblue')
        
        for arg_key in args_provided:
            # Extract input value into separate yml file
            # Replace it here with a new variable name
            arg_val = steps[i][key]['in'][arg_key]
            var_name = f'{step_name}/{arg_key}'
            in_name = f'{step_name}_input___{arg_key}'  # Use triple underscore for namespacing so we can split later
            if arg_val[0] == '$':
                #print('arg_key, arg_val', arg_key, arg_val)
                if not vars_workflow_input.get(arg_val[1:]):
                    inputs_workflow.update({in_name: arg_val[1:]})
                    inputs_file_workflow.update({in_name: arg_val[1:]})
                    steps[i][key]['in'][arg_key] = in_name
                    vars_workflow_input.update({arg_val[1:]: var_name})
                    if args.graph_show_inputs:
                        dot.node(var_name, label=arg_key, shape='box', style='rounded, filled', fillcolor='lightgreen')
                        dot.edge(var_name, step_name)
                else:
                    steps[i][key]['in'][arg_key] = vars_workflow_input[arg_val[1:]]
                    dot.edge(vars_workflow_input[arg_val[1:]].split('/')[0], step_name)
            else:
                inputs_workflow.update({in_name: arg_val})
                inputs_file_workflow.update({in_name: arg_val})
                steps[i][key]['in'][arg_key] = in_name
                if args.graph_show_inputs:
                    dot.node(var_name, label=arg_key, shape='box', style='rounded, filled', fillcolor='lightgreen')
                    dot.edge(var_name, step_name)
        
        for arg_key in args_required:
            #print('arg_key', arg_key)
            if arg_key in args_provided:
                continue  # We already convered this case above.
            match = False
            # TODO: Figure out something better than replace
            arg_key_no_namespace = arg_key.split('___')[-1]
            arg_key_noinput = arg_key_no_namespace.replace('input_', '').replace('solute_', '').replace('energy_', 'edr_').replace('structure_', 'tpr_').replace('traj_', 'trr_')
            # gmx_trjconv_str assumes gro, but other tools assume tpr; For now,
            # instead of figuring out how to replace('structure_', 'gro_') for
            # gmx_trjconv_str only, just require the user to specify filename.
            for j in range(0, i)[::-1]:  # Reverse order!
                out_keys = list(tools[j]['outputs'])[::-1] # Reverse order!
                for out_key in out_keys:
                    out_key_no_namespace = out_key.split('___')[-1]
                    if arg_key_noinput == out_key_no_namespace.replace('output_', ''):
                        #print('match!', j)  # We found a match!
                        # Generate a new namespace for out_key using the step number and add to inputs
                        step_name_j = f'step_{j + 1}_{steps_keys[j]}'
                        arg_val = f'{step_name_j}/{out_key}'
                        arg_keyval = {arg_key: arg_val}
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

                        # Determine which tail node to use for the new edge
                        edge_node1 = step_name_j
                        if tools[j]['class'] == 'Workflow' and steps_keys[j] in subkeys: # i.e. if we performed recursion
                            if args.graph_show_outputs:
                                if args.graph_inline_subgraphs:
                                    edge_node1 = out_key.split('___')[0] # out_key
                            else:
                                edge_node1 = out_key.split('___')[0]

                        # We also need to keep track of the 'internal' output variables
                        if tools[j]['class'] == 'Workflow':
                            vars_workflow_output_internal.append(out_key)
                        else:
                            vars_workflow_output_internal.append(arg_val)
                        
                        # Add an edge between two steps in the graph
                        if args.graph_label_edges:
                            if tools[j]['class'] == 'Workflow':
                                dot.edge(edge_node1, step_name, label=out_key_no_namespace)
                            else:
                                dot.edge(edge_node1, step_name, label=out_key)
                        else:
                            dot.edge(edge_node1, step_name)
                # Break out two levels, i.e. continue onto the next iteration of the outermost loop.
                        match = True
                        break
                if match:
                    break
            if not match:
                # This just means we need to defer to the parent workflow.
                # There will actually be an error only if no parent supplies an input value.
                if is_root:
                    print('Error! No match found for input', i + 1, key, arg_key)

                in_name = f'{step_name}_input___{arg_key}'  # Use triple underscore for namespacing so we can split later
                arg_keyval = {arg_key: in_name}
                # Add an input name to this subworkflow. Add a namespace prefix in the parent workflow.
                inputs_workflow.update({in_name: arg_key})
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
        out_keys = list(tools[i]['outputs'])
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

    # Add the provided inputs of each step to the workflow inputs
    def filetype_predicate(key):
        # TODO: Improve File type detection heuristics
        keys = key.split('___')
        return 'path' in keys[-1] and (not 'xvg' in keys[-1]) and ('grompp' in keys[0] and not 'tpr' in keys[-1]) and (not 'pdb' in keys[0])
    inputs_workflow_types = dict([(key, 'File' if filetype_predicate(key) else 'string') for key in inputs_workflow])
    yaml_tree.update({'inputs': inputs_workflow_types})

    # Add the outputs of each step to the workflow outputs
    vars_workflow_output_internal = list(set(vars_workflow_output_internal))  # Get uniques
    outputs_workflow = {}
    for i, key in enumerate(steps_keys):
        step_name = f'step_{i + 1}_{steps_keys[i]}'
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
                case1 = (tools[i]['class'] == 'Workflow') and (not out_key in [var.replace('/', '___') for var in vars_workflow_output_internal])
                case2 = (tools[i]['class'] == 'CommandLineTool') and (not out_var in vars_workflow_output_internal)
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
        step_name = f'step_{i + 1}_{steps_keys[i]}'
        #steps[i] = {step_name: steps[i][key]}
        steps_dict.update({step_name: steps[i][key]})
    yaml_tree.update({'steps': steps_dict})
    
    # Dump the workflow inputs to a separate yml file.
    for key, val in inputs_file_workflow.items():
        if filetype_predicate(key):
            inputs_file_workflow.update({key: {'class': 'File', 'path': val, 'format': 'https://edamontology.org/format_2033'}})
    dump_options = {'line_break': '\n', 'indent': 2}
    yaml_content = yaml.dump(inputs_file_workflow, sort_keys=False, **dump_options)
    with open(f'{yaml_stem}_inputs.yml', 'w') as inp:
        inp.write(yaml_content)

    # Dump the expanded yaml file to disk.
    # Use sort_keys=False to preserve the order of the steps.
    dump_options = {'line_break': '\n', 'indent': 2}
    yaml_content = yaml.dump(yaml_tree, sort_keys=False, **dump_options)
    with open(f'{yaml_stem}.cwl', 'w') as w:
        w.write('#!/usr/bin/env cwl-runner\n')
        w.write(''.join(yaml_content))
    
    print('Finished expanding. Validating...')
    
    cmd = ['cwltool', '--validate', f'{yaml_stem}.cwl']
    sub.run(cmd)
    
    return (yaml_tree, inputs_file_workflow, vars_workflow_output_internal, is_root, dot)


def main():
    parser = argparse.ArgumentParser(prog='main', description='Convert a high-level yaml workflow file to CWL.')
    parser.add_argument('--yaml', type=str, required=True,
                        help='Yaml workflow file')
    parser.add_argument('--clt_dir', type=str, required=False, default='biobb',
                        help='Directory which contains the CWL CommandLineTools')
    parser.add_argument('--cwl_output_intermediate_files', type=bool, required=False, default=False,
                        help='Enable output files which are used between steps (for debugging).')
    parser.add_argument('--graph_label_edges', type=bool, required=False, default=False,
                        help='Label the graph edges with the name of the intermediate input/output.')
    parser.add_argument('--graph_show_inputs', type=bool, required=False, default=False,
                        help='Add nodes to the graph representing the workflow inputs.')
    parser.add_argument('--graph_show_outputs', type=bool, required=False, default=False,
                        help='Add nodes to the graph representing the workflow outputs.')
    parser.add_argument('--graph_inline_subgraphs', type=bool, required=False, default=False,
                        help='Controls whether subgraphs are displayed separately or positioned within the main graph.')
    args = parser.parse_args()

    # Load ALL of the tools.
    tools_cwl = {}
    pattern_cwl = str(Path(args.clt_dir) / '**/*.cwl')
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

    dot_in = graphviz.Digraph(name=args.yaml)
    #dot_in.attr(rankdir='LR')
    workflow_data = expand_workflow(args, dot_in, tools_cwl, True, args.yaml)
    dot_out = workflow_data[-1]
    # Render the GraphViz diagram
    dot_out.render()
    #dot_out.view() # viewing does not work on headless machines (and requires xdg-utils)


if __name__ == '__main__':
    main()
