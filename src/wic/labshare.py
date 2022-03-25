import requests
import yaml


def upload_plugin(compute_url, tool, stem):
        # Convert the compiled yaml file to json for labshare Compute.
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


def print_plugins(compute_url):
    r = requests.get(compute_url + '/compute/plugins/')
    for j in r.json():
        print(j)
    print(len(r.json()))