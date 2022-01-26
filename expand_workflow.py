
import argparse
import copy
import glob
from pathlib import Path
import subprocess as sub

import yaml

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='main', description='Convert a high-level yaml workflow file to CWL.')
    parser.add_argument('--yaml', type=str, required=True,
                        help='Yaml workflow file')
    parser.add_argument('--clt_dir', type=str, required=False, default='biobb',
                        help='Directory which contains the CWL CommandLineTools')
    args = parser.parse_args()

    # Load the high-level yaml workflow file.
    yaml_stem = Path(args.yaml).stem
    with open(Path(args.yaml), 'r') as y:
        yaml_tree = yaml.load(y.read(), Loader=yaml.SafeLoader)

    #for x in yaml_tree:
    #    print(x)
    steps = yaml_tree['steps']
    #for step in steps:
    #  print(step)
    
    steps_keys = []
    for step in steps:
        steps_keys += list(step)
    #print(steps_keys)

    # Load ALL of the tools.
    tools_all = {}
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
            tools_all[stem] = (cwl_path, tool)
            #print(tool)
        except yaml.scanner.ScannerError as se:
            pass
            # There are two cwl files that throw this error, but they are both legacy, so...
            #print(cwl_path)
            #print(se)

    # Restrict to the subset of tools that we are actually using.
    tools = [tools_all[key][1] for key in steps_keys]

    # Collect workflow input parameters
    inputs_workflow = {}

    for i, key in enumerate(steps_keys):
        # Add run tag
        if steps[i][key]:
            if not 'run' in steps[i][key]:
                steps[i][key].update({'run': tools_all[key][0]})
        else:
            steps[i] = {key: {'run': tools_all[key][0]}}

        # Generate intermediate file names between steps.
        # The inputs for a given step need to come from either
        # 1. the input.yaml file for the overall workflow (extract into separate yml file) or
        # 2. the outputs from the previous steps (most recent first).
        # If there isn't an exact match, remove input_* and output_* from the
        # current step and previous steps, respectively, and then check again.
        # If there still isn't an exact match, explicit renaming may be required.

        in_tool = tools[i]['inputs']
        #print(list(in_tool.keys()))
        args_required = [arg for arg in in_tool if not (in_tool[arg].get('default') or in_tool[arg]['type'][-1] == '?')]
        # Note: Some config tags are not required in the cwl files, but are in
        # fact required in the python source code! See check_mandatory_property
        # (Solution: refactor all required arguments out of config and list
        # them as explicit inputs in the cwl files, then modify the python
        # files accordingly.)
        #print(args_required)

        args_provided = []
        if steps[i][key] and 'in' in steps[i][key]:
            args_provided = list(steps[i][key]['in'])
        #print(args_provided)
        
        step_name = f'step_{i + 1}_{steps_keys[i]}'
        
        for arg_key in args_provided:
            # Extract input value into separate yml file
            # Replace it here with a new variable name
            var_name = f'{step_name}_{arg_key}'
            inputs_workflow.update({var_name: steps[i][key]['in'][arg_key]})
            steps[i][key]['in'][arg_key] = var_name
        
        for arg_key in args_required:
            #print('arg_key', arg_key)
            if arg_key in args_provided:
                continue  # We already convered this case above.
            match = False
            # TODO: Figure out something better than replace
            arg_key_noinput = arg_key.replace('input_', '').replace('solute_', '').replace('energy_', 'edr_').replace('structure_', 'tpr_').replace('traj_', 'trr_')
            # gmx_trjconv_str assumes gro, but other tools assume tpr; For now,
            # instead of figuring out how to replace('structure_', 'gro_') for
            # gmx_trjconv_str only, just require the user to specify filename.
            for j in range(0, i)[::-1]:  # Reverse order!
                out_keys = list(tools[j]['outputs'])
                for out_key in out_keys:
                    if out_key.replace('output_', '') == arg_key_noinput:
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
                # Break out two levels, i.e. continue onto the next iteration of the outermost loop.
                        match = True
                        break
                if match:
                    break
            if not match:
                print('Error! No match found for input arg', arg_key)
        
        # Add CommandLineTool outputs tags to workflow out tags.
        # Note: Add all output tags for now, but depending on config options,
        # not all output files will be genenerated. This may cause an error.
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
    inputs_workflow_types = dict([(key, 'File' if 'path' in key and not 'xvg' in key and not 'step_1' in key and not 'step_8' in key else 'string') for key in inputs_workflow])
    yaml_tree.update({'inputs': inputs_workflow_types})

    # Add the outputs of each step to the workflow outputs
    outputs_workflow = {}
    for i, key in enumerate(steps_keys):
        step_name = f'step_{i + 1}_{steps_keys[i]}'
        out_keys = steps[i][key]['out']
        for out_key in out_keys:
            # Exclude certain output files as per the comment above.
            if 'dhdl' in out_key or 'xtc' in out_key:
                continue
            out_var = f'{step_name}/{out_key}'
            out_name = f'{step_name}_{out_key}'
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
    for key, val in inputs_workflow.items():
        # TODO: Improve File type detection heuristics
        if 'path' in key and not 'xvg' in key and not 'step_1' in key and not 'step_8' in key :
            inputs_workflow.update({key: {'class': 'File', 'path': val, 'format': 'https://edamontology.org/format_2033'}})
    dump_options = {'line_break': '\n', 'indent': 2}
    yaml_content = yaml.dump(inputs_workflow, sort_keys=False, **dump_options)
    with open(f'{yaml_stem}_inputs.yml', 'w') as inp:
        inp.write(yaml_content)

    # Dump the expanded yaml file to disk.
    # Use sort_keys=False to preserve the order of the steps.
    dump_options = {'line_break': '\n', 'indent': 2}
    yaml_content = yaml.dump(yaml_tree, sort_keys=False, **dump_options)
    with open(f'{yaml_stem}_expanded.cwl', 'w') as w:
        w.write('#!/usr/bin/env cwl-runner\n')
        w.write('cwlVersion: v1.0\n')
        w.write('class: Workflow\n')
        w.write(''.join(yaml_content))
    
    print('Finished expanding. Validating...')
    
    cmd = ['cwltool', '--validate', f'{yaml_stem}_expanded.cwl']
    sub.run(cmd)