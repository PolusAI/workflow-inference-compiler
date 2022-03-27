from typing import Dict, Any

import requests
import yaml

from .wic_types import KV, Cwl

def upload_plugin(compute_url: str, tool: Cwl, stem: str) -> int:
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

        plugin_id: str = r_json['id']
        compute_plugin['id'] = plugin_id
        compute_plugin.update({'id': plugin_id}) # Necessary ?
        # Save the plugin ids so we can delete them the next time we enable args.cwl_run_slurm
        with open('plugin_ids', 'a') as f:
            f.write(plugin_id + '\n')
        return int(plugin_id)  # hash


def print_plugins(compute_url: str) -> None:
    r = requests.get(compute_url + '/compute/plugins/')
    for j in r.json():
        print(j)
    print(len(r.json()))