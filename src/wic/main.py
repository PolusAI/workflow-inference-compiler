import glob
from pathlib import Path
import requests
import subprocess as sub

import graphviz
import yaml

from . import cli
from . import compiler
from . import utils

# Use white for dark backgrounds, black for light backgrounds
font_edge_color = 'white'


def main():
    args = cli.parser.parse_args()

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
        steps_inlined = utils.inline_sub_steps(yaml_path, tools_cwl, yml_paths)
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
        workflow_data = compiler.compile_workflow(args, [], [subgraph], vars_dollar_defs, tools_cwl, True, yaml_path, yml_paths)
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