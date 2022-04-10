import argparse
from pathlib import Path
from typing import Dict, List

import requests
import yaml

from .wic_types import KV, Cwl, NodeData, RoseTree, Tools
from . import utils

def upload_plugin(compute_url: str, tool: Cwl, stem: str) -> str:
    """Uploads CWL CommandLineTools to Polus Compute

    Args:
        compute_url (str): The url to the Compute API
        tool (Cwl): The CWL CommandLineTool
        stem (str): The name of the CWL CommandLineTool

    Raises:
        Exception: If the upload failed for any reason

    Returns:
        str: The unique id of the plugin
    """
    # Convert the compiled yaml file to json for labshare Compute.
    # First remove $ in $namespaces and $schemas (anywhere else?)
    tool_str = str(yaml.dump(tool, sort_keys=False, line_break='\n', indent=2))
    tool_str_no_dollar = tool_str.replace('$namespaces', 'namespaces').replace('$schemas', 'schemas')
    tool_no_dollar: Cwl = yaml.safe_load(tool_str_no_dollar)  # This effectively copies tool
    compute_plugin: KV = {
        # Add unique 'id' below
        'cwlScript': tool_no_dollar
    }

    # Use http POST request to upload a primitive CommandLineTool / define a plugin and get its id hash.
    response = requests.post(compute_url + '/compute/plugins', json = compute_plugin)
    r_json = response.json()
    if 'id' not in r_json:
        print('post response')
        print(r_json)
        raise Exception(f'Error! Labshare plugin upload failed for {stem}.')

    plugin_id: str = r_json['id'] # hash
    compute_plugin['id'] = plugin_id
    compute_plugin.update({'id': plugin_id}) # Necessary ?
    # Save the plugin ids so we can delete them the next time we enable args.cwl_run_slurm
    with open('plugin_ids', 'a') as f:
        f.write(f'{plugin_id}\n')
    return plugin_id


def print_plugins(compute_url: str) -> None:
    """prints information on all currently available Compute plugins

    Args:
        compute_url (str): The url to the Compute API
    """
    r = requests.get(compute_url + '/compute/plugins/')
    for j in r.json():
        print(f"id {j.get('id')} class {j.get('class')} name {j.get('name')}")
        #print(j)
    print(len(r.json()))


def upload_all(rose_tree: RoseTree, tools: Tools, args: argparse.Namespace, is_root: bool) -> str:
    """Uploads all Plugins, Pipelines, and the root Workflow to the Compute platform

    Args:
        rose_tree (RoseTree): The data associated with compiled subworkflows
        tools (Tools): The CWL CommandLineTool definitions found using get_tools_cwl()
        args (argparse.Namespace): The command line arguments
        is_root (bool): True if this is the root workflow

    Raises:
        Exception: If any of the uploads fails for any reason

    Returns:
        str: The unique id of the workflow
    """
    sub_node_data: NodeData = rose_tree.data
    yaml_stem = sub_node_data.name
    cwl_tree = sub_node_data.compiled_cwl
    yaml_inputs = sub_node_data.workflow_inputs_file
    
    sub_rose_trees: Dict[str, RoseTree] = dict([(r.data.name, r) for r in rose_tree.sub_trees])
    #print(list(sub_rose_trees))

    steps = cwl_tree['steps']

    # Get the dictionary key (i.e. the name) of each step.
    steps_keys: List[str] = []
    for step in steps:
        step_key = utils.parse_step_name_str(step)[-1]
        steps_keys.append(step_key)
    #print(steps_keys)

    #subkeys = [key for key in steps_keys if key not in tools]

    # Convert the compiled yaml file to json for labshare Compute.
    # Replace 'run' with plugin:id
    import copy
    cwl_tree_run = copy.deepcopy(cwl_tree)
    for i, step_key in enumerate(steps_keys):
        stem = Path(step_key).stem
        tool_i = tools[stem].cwl
        step_name_i = utils.step_name_str(yaml_stem, i, step_key)

        #if step_key in subkeys: # and not is_root, but the former implies the latter
            #plugin_id = upload_plugin(args.compute_url, cwl_tree_run, yaml_stem)
        if stem in sub_rose_trees:
            subworkflow_id = upload_all(sub_rose_trees[stem], tools, args, False)
            run_val = f'pipeline:{subworkflow_id}'

            # Save the pipeline ids so we can delete them the next time we enable args.cwl_run_slurm
            with open('pipeline_ids', 'a') as f:
                f.write(f'{subworkflow_id}\n')
        else:
            # i.e. If this is either a primitive CommandLineTool and/or
            # a 'primitive' Workflow that we did NOT recursively generate.
            plugin_id = upload_plugin(args.compute_url, tool_i, stem)
            run_val = f'plugin:{plugin_id}'
        cwl_tree_run['steps'][step_name_i]['run'] = run_val

    workflow_id: str = ''
    if is_root:
        compute_workflow = {
            "name": yaml_stem,
            "driver": "slurm",
            "cwlJobInputs": yaml_inputs,
            **cwl_tree_run
        }
        # Use http POST request to upload a complete Workflow (w/ inputs) and get its id hash.
        response = requests.post(args.compute_url + '/compute/workflows', json = compute_workflow)
        r_json = response.json()
        print('post response')
        j = r_json
        print(f"id {j.get('id')} class {j.get('class')} name {j.get('name')}")
        if 'id' not in r_json:
            print(r_json)
            raise Exception(f'Error! Labshare workflow upload failed for {yaml_stem}.')
        workflow_id = r_json['id'] # hash
    else:
        #  "owner": "string",
        #  "additionalProp1": {}
        # TODO: Check this.
        compute_pipeline = {
            "name": yaml_stem,
            **cwl_tree_run
        }
        # Need to add owner and/or additionalProp1 ?
        # Need to remove headers and/or requirements? i.e.
        #yaml_tree['cwlVersion'] = 'v1.0'
        #yaml_tree['class'] = 'Workflow'
        #yaml_tree['requirements'] = subworkreqdict
    
        # Use http POST request to upload a subworkflow / "pipeline" (no inputs) and get its id hash.
        response = requests.post(args.compute_url + '/compute/pipelines', json = compute_pipeline)
        r_json = response.json()
        print('post response')
        j = r_json
        print(f"id {j.get('id')} class {j.get('class')} name {j.get('name')}")
        if 'id' not in r_json:
            print(r_json)
            raise Exception(f'Error! Labshare workflow upload failed for {yaml_stem}.')
        workflow_id = r_json['id'] # hash
    if is_root:
        print_plugins(args.compute_url)

    return workflow_id