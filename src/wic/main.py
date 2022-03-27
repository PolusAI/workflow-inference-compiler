import glob
from pathlib import Path
import requests
import subprocess as sub
from typing import Dict

import graphviz
import yaml

from . import auto_gen_header
from . import cli
from . import compiler
from . import labshare
from . import utils
from .wic_types import Cwl, Yaml, Tools

# Use white for dark backgrounds, black for light backgrounds
font_edge_color = 'white'


def get_tools_cwl(cwl_dir: Path) -> Tools:
    # Load ALL of the tools.
    tools_cwl: Tools = {}
    pattern_cwl = str(cwl_dir / '**/*.cwl')
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

    for cwl_path_str in cwl_paths_sorted:
        #print(cwl_path)
        try:
            with open(cwl_path_str, 'r') as f:
              tool: Cwl = yaml.safe_load(f.read())
            stem = Path(cwl_path_str).stem
            #print(stem)
            # Add / overwrite stdout and stderr
            tool.update({'stdout': f'{stem}.out'})
            tool.update({'stderr': f'{stem}.err'})
            tools_cwl[stem] = (cwl_path_str, tool)
            #print(tool)
        except yaml.scanner.ScannerError as se:
            pass
            # There are two cwl files that throw this error, but they are both legacy, so...
            #print(cwl_path)
            #print(se)
    
    return tools_cwl


def get_yml_paths(cwl_dir: Path) -> Dict[str, Path]:
    # Glob all of the yml files too, so we don't have to deal with relative paths.
    pattern_yml = str(cwl_dir / '**/*.yml')
    yml_paths_sorted = sorted(glob.glob(pattern_yml, recursive=True), key=len, reverse=True)
    yml_paths = {}
    for yml_path_str in yml_paths_sorted:
        yml_path = Path(yml_path_str)
        yml_paths[yml_path.stem] = yml_path

    return yml_paths


def main() -> None:
    args = cli.parser.parse_args()

    tools_cwl = get_tools_cwl(Path(args.cwl_dir))
    yml_paths = get_yml_paths(Path(args.cwl_dir))

    yaml_path = args.yaml
    if args.cwl_inline_subworkflows:
        steps_inlined = utils.inline_sub_steps(yaml_path, tools_cwl, yml_paths)
        with open(Path(yaml_path), 'r') as y:
            yaml_tree: Yaml = yaml.safe_load(y.read())
        yaml_tree['steps'] = steps_inlined
        yaml_content = yaml.dump(yaml_tree, sort_keys=False, line_break='\n', indent=2)
        yaml_path = Path(args.yaml).stem + '_inline.yml'
        with open(yaml_path, 'w') as y:
            y.write(auto_gen_header)
            y.write(yaml_content)

    rootgraph = graphviz.Digraph(name=yaml_path)
    rootgraph.attr(newrank='True') # See graphviz layout comment above.
    rootgraph.attr(bgcolor="transparent") # Useful for making slides
    rootgraph.attr(fontcolor=font_edge_color)
    #rootgraph.attr(rankdir='LR') # When --graph_inline_depth 1, this usually looks better.
    with rootgraph.subgraph(name=f'cluster_{yaml_path}') as subgraph:
        subgraph.attr(label=yaml_path)
        subgraph.attr(color='lightblue')  # color of cluster subgraph outline
        workflow_data = compiler.compile_workflow(args, [], [subgraph], {}, tools_cwl, True, yaml_path, yml_paths)
        recursive_data = workflow_data[0]

    utils.write_to_disk(recursive_data)

    if args.cwl_run_slurm:
        # Delete plugins previously uploaded to labshare.
        if args.cwl_run_slurm and Path('plugin_ids').exists():
            with open('plugin_ids', 'r') as f:
                ids = f.read().splitlines()
            for id in ids:
                response = requests.delete(args.compute_url + '/compute/plugins/' + id)
            sub.run(['rm', 'plugin_ids'])
        # Delete pipelines previously uploaded to labshare.
        if args.cwl_run_slurm and Path('pipeline_ids').exists():
            with open('pipeline_ids', 'r') as f:
                ids = f.read().splitlines()
            for id in ids:
                response = requests.delete(args.compute_url + '/compute/pipelines/' + id)
            sub.run(['rm', 'pipeline_ids'])
        labshare.upload_all(recursive_data, tools_cwl, args, True)

    # Render the GraphViz diagram
    rootgraph.render(format='png') # Default pdf. See https://graphviz.org/docs/outputs/
    #rootgraph.view() # viewing does not work on headless machines (and requires xdg-utils)
    
    if args.cwl_run_local:
        yaml_stem = Path(args.yaml).stem
        yaml_stem = yaml_stem + '_inline' if args.cwl_inline_subworkflows else yaml_stem
        print(f'Running {yaml_stem}.cwl ...')
        cmd = ['cwltool', '--cachedir', 'cachedir','--outdir', 'outdir', f'{yaml_stem}.cwl', f'{yaml_stem}_inputs.yml']
        sub.run(cmd)


if __name__ == '__main__':
    main()