import glob
import subprocess as sub
import sys
import os
from pathlib import Path
from typing import Dict

import graphviz
import networkx as nx
import yaml

from . import ast, cli, compiler, inference, labshare, utils
from .schemas import wic_schema
from .wic_types import Cwl, GraphData, GraphReps, StepId, Tool, Tools, Yaml, YamlTree


def get_tools_cwl(cwl_dirs_file: Path) -> Tools:
    """Uses glob() to find all of the CWL CommandLineTool definition files within any subdirectory of cwl_dir

    Args:
        cwl_dirs_file (Path): The subdirectories in which to search for CWL CommandLineTools

    Returns:
        Tools: The CWL CommandLineTool definitions found using glob()
    """
    utils.copy_config_files()
    # Load ALL of the tools.
    tools_cwl: Tools = {}
    cwl_dirs = utils.read_lines_pairs(cwl_dirs_file)
    for plugin_ns, cwl_dir in cwl_dirs:
        # "PurePath.relative_to() requires self to be the subpath of the argument, but os.path.relpath() does not."
        # See https://docs.python.org/3/library/pathlib.html#id4 and
        # See https://stackoverflow.com/questions/67452690/pathlib-path-relative-to-vs-os-path-relpath
        cwl_dir_rel = os.path.relpath(cwl_dir) # w.r.t. current working directory
        pattern_cwl = str(Path(cwl_dir_rel) / '**/*.cwl')
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
        Path('autogenerated/schemas/tools/').mkdir(parents=True, exist_ok=True)
        if len(cwl_paths_sorted) == 0:
            print(f'Warning! No cwl files found in {cwl_dir_rel}.\nCheck {cwl_dirs_file.absolute()}')
            print('This almost certainly means you are not in the correct directory.')
        for cwl_path_str in cwl_paths_sorted:
            #print(cwl_path)
            try:
                with open(cwl_path_str, mode='r', encoding='utf-8') as f:
                    tool: Cwl = yaml.safe_load(f.read())
                stem = Path(cwl_path_str).stem
                # print(stem)
                # Add / overwrite stdout and stderr
                tool.update({'stdout': f'{stem}.out'})
                tool.update({'stderr': f'{stem}.err'})
                step_id = StepId(stem, plugin_ns)
                tools_cwl[step_id] = Tool(cwl_path_str, tool)
                #print(tool)
            except yaml.scanner.ScannerError as s_e:
                pass
                # There are two cwl files that throw this error, but they are both legacy, so...
                #print(cwl_path)
                #print(s_e)
            #utils.make_tool_dag(stem, (cwl_path_str, tool))
    return tools_cwl


def get_yml_paths(yml_dirs_file: Path) -> Dict[str, Dict[str, Path]]:
    """Uses glob() to recursively find all of the yml workflow definition files
    within any subdirectory of each yml_dir in yml_dirs_file.
    NOTE: This function assumes all yml files found are workflow definition files,
    so do not mix regular yml files and workflow files in the same root directory.
    Moreover, each yml_dir should be disjoint; do not use both '.' and './subdir'!

    Args:
        yml_dirs_file (Path): The subdirectories in which to search for yml files

    Returns:
        Dict[str, Dict[str, Path]]: A dict containing the filepath stem and filepath of each yml file
    """
    utils.copy_config_files()
    yml_dirs = utils.read_lines_pairs(yml_dirs_file)
    # Glob all of the yml files too, so we don't have to deal with relative paths.
    yml_paths_all: Dict[str, Dict[str, Path]] = {}
    for yml_namespace, yml_dir in yml_dirs:
        # "PurePath.relative_to() requires self to be the subpath of the argument, but os.path.relpath() does not."
        # See https://docs.python.org/3/library/pathlib.html#id4 and
        # See https://stackoverflow.com/questions/67452690/pathlib-path-relative-to-vs-os-path-relpath
        yml_dir_rel = os.path.relpath(yml_dir) # w.r.t. current working directory
        pattern_yml = str(Path(yml_dir_rel) / '**/*.yml')
        yml_paths_sorted = sorted(glob.glob(pattern_yml, recursive=True), key=len, reverse=True)
        if len(yml_paths_sorted) == 0:
            print(f'Warning! No yml files found in {yml_dir_rel}.\nCheck {yml_dirs_file.absolute()}')
            print('This almost certainly means you are not in the correct directory.')
        yml_paths = {}
        for yml_path_str in yml_paths_sorted:
            # Exclude our autogenerated inputs files
            if '_inputs' not in yml_path_str:
                yml_path = Path(yml_path_str)
                yml_paths[yml_path.stem] = yml_path
        # Check for existing entry (so we can split a single
        # namespace across multiple lines in yml_dirs.txt)
        ns_dict = yml_paths_all.get(yml_namespace, {})
        yml_paths_all[yml_namespace] = {**ns_dict, **yml_paths}

    return yml_paths_all


def main() -> None:
    """See docs/userguide.md"""
    args = cli.parser.parse_args()

    tools_cwl = get_tools_cwl(args.cwl_dirs_file)
    utils.make_plugins_dag(tools_cwl, args.graph_dark_theme)
    yml_paths = get_yml_paths(args.yml_dirs_file)

    # Perform initialization via mutating global variables (This is not ideal)
    compiler.inference_rules = dict(utils.read_lines_pairs(Path('inference_rules.txt')))
    inference.renaming_conventions = utils.read_lines_pairs(Path('renaming_conventions.txt'))

    # Generate schemas for validation and vscode IntelliSense code completion
    yaml_stems = utils.flatten([list(p) for p in yml_paths.values()])
    validator = wic_schema.get_validator(tools_cwl, yaml_stems, write_to_disk=True)
    if args.generate_schemas_only:
        print('Finished generating schemas. Exiting.')
        sys.exit(0)

    yaml_path = args.yaml

    # Load the high-level yaml root workflow file.
    with open(yaml_path, mode='r', encoding='utf-8') as y:
        root_yaml_tree: Yaml = yaml.safe_load(y.read())
    Path('autogenerated/').mkdir(parents=True, exist_ok=True)
    wic = {'wic': root_yaml_tree.get('wic', {})}
    plugin_ns = wic['wic'].get('namespace', 'global')
    step_id = StepId(yaml_path, plugin_ns)
    y_t = YamlTree(step_id, root_yaml_tree)
    yaml_tree_raw = ast.read_ast_from_disk(y_t, yml_paths, tools_cwl, validator)
    # Write the combined workflow (with all subworkflows as children) to disk.
    with open(f'autogenerated/{Path(yaml_path).stem}_tree_raw.yml', mode='w', encoding='utf-8') as f:
        f.write(yaml.dump(yaml_tree_raw.yml))
    yaml_tree = ast.merge_yml_trees(yaml_tree_raw, {}, tools_cwl)
    with open(f'autogenerated/{Path(yaml_path).stem}_tree_merged.yml', mode='w', encoding='utf-8') as f:
        f.write(yaml.dump(yaml_tree.yml))

    if args.cwl_inline_subworkflows:
        while True:
            # Inlineing changes the namespaces, so we have to get new namespaces after each inlineing operation.
            namespaces_list = ast.get_inlineable_subworkflows(yaml_tree, tools_cwl, False, [])
            if namespaces_list == []:
                break

            #print('inlineing', namespaces_list[0])
            yaml_tree = ast.inline_subworkflow(yaml_tree, tools_cwl, namespaces_list[0])

        with open(f'autogenerated/{Path(yaml_path).stem}_tree_merged_inlined.yml', mode='w', encoding='utf-8') as f:
            f.write(yaml.dump(yaml_tree.yml))

    rootgraph = graphviz.Digraph(name=yaml_path)
    rootgraph.attr(newrank='True') # See graphviz layout comment above.
    rootgraph.attr(bgcolor="transparent") # Useful for making slides
    font_edge_color = 'black' if args.graph_dark_theme else 'white'
    rootgraph.attr(fontcolor=font_edge_color)

    # This can be used to visually 'inline' all subworkflows (but NOT the CWL).
    # rootgraph.attr(style='invis')
    # Note that since invisible objects still affect the graphviz layout (by design),
    # this can be used to control the layout of the individual nodes, even if
    # you don't necessarily want subworkflows.

    #rootgraph.attr(rankdir='LR') # When --graph_inline_depth 1, this usually looks better.
    with rootgraph.subgraph(name=f'cluster_{yaml_path}') as subgraph_gv:
        # get the label (if any) from the workflow
        step_i_wic_graphviz = yaml_tree.yml.get('wic', {}).get('graphviz', {})
        label = step_i_wic_graphviz.get('label', yaml_path)
        subgraph_gv.attr(label=label)
        subgraph_gv.attr(color='lightblue')  # color of cluster subgraph outline
        subgraph_nx = nx.DiGraph()
        graphdata = GraphData(yaml_path)
        subgraph = GraphReps(subgraph_gv, subgraph_nx, graphdata)
        compiler_info = compiler.compile_workflow(yaml_tree, args, [], [subgraph], {}, {},
                                                  tools_cwl, True, relative_run_path=True, testing=False)
        rose_tree = compiler_info.rose

    utils.write_to_disk(rose_tree, Path('autogenerated/'), relative_run_path=True)

    if args.run_compute:
        # Inline compiled CWL if necessary, i.e. inline across scattering boundaries.
        # NOTE: Since we need to distribute scattering operations across all dependencies,
        # and due to inference, this cannot be done before compilation.
        rose_tree = ast.inline_subworkflow_cwl(rose_tree)
        utils.write_to_disk(rose_tree, Path('autogenerated/'), relative_run_path=True)
        labshare.upload_all(rose_tree, tools_cwl, args, True)

    # Render the GraphViz diagram
    rootgraph.render(format='png') # Default pdf. See https://graphviz.org/docs/outputs/
    yaml_stem = Path(args.yaml).stem
    #cmd = f'cwltool --print-dot autogenerated/{yaml_stem}.cwl | dot -Tsvg > autogenerated/{yaml_stem}.svg'
    #sub.run(cmd, shell=True, check=False)
    #rootgraph.view() # viewing does not work on headless machines (and requires xdg-utils)

    if args.run_local:
        yaml_inputs = rose_tree.data.workflow_inputs_file
        stage_input_files(yaml_inputs, Path(args.yaml).parent.absolute())

        yaml_stem = yaml_stem + '_inline' if args.cwl_inline_subworkflows else yaml_stem
        if args.cwl_runner == 'cwltool':
            parallel = ['--parallel'] if args.parallel else []
            # NOTE: --parallel is required for real-time analysis / real-time plots,
            # but it seems to cause hanging with Docker for Mac. The hanging seems
            # to be worse when using parallel scattering.
            quiet = ['--quiet'] if args.quiet else []
            # NOTE: Using --leave-outputs to disable --outdir
            # See https://github.com/dnanexus/dx-cwl/issues/20
            # --outdir has one or more bugs which will cause workflows to fail!!!
            cmd = ['cwltool'] + parallel + quiet + ['--leave-tmpdir', '--leave-outputs',
                    '--provenance', 'provenance', '--cachedir', args.cachedir,
                f'autogenerated/{yaml_stem}.cwl', f'autogenerated/{yaml_stem}_inputs.yml']
            # '--js-console' "Running with support for javascript console in expressions (DO NOT USE IN PRODUCTION)"
            # TODO: Consider using the undocumented flag --fast-parser for known-good workflows,
            # which was recently added in the 3.1.20220913185150 release of cwltool.
        if args.cwl_runner == 'toil-cwl-runner':
            # NOTE: toil-cwl-runner always runs in parallel
            cmd = ['toil-cwl-runner', '--provenance', 'provenance', '--outdir', 'outdir_toil',
                '--jobStore', 'file:./jobStore', # NOTE: This is the equivalent of --cachedir
                # TODO: Check --clean, --cleanWorkDir, --restart
                f'autogenerated/{yaml_stem}.cwl', f'autogenerated/{yaml_stem}_inputs.yml']
        print('Running ' + ' '.join(cmd))
        proc = sub.run(cmd, check=False)
        if proc.returncode == 0:
            print('Success! Output files should be in outdir/')
        else:
            print('Failure! Please scroll up and find the FIRST error message.')
            print('(You may have to scroll up A LOT.)')

        # Finally, since there is an output file copying bug in cwltool,
        # we need to copy the output files manually. See comment above.
        output_json_file = Path('provenance/workflow/primary-output.json')
        if output_json_file.exists():
            files = utils.parse_provenance_output_files(output_json_file)

            dests = set()
            for location, namespaced_output_name, basename in files:
                yaml_stem_init, shortened = utils.shorten_namespaced_output_name(namespaced_output_name)
                parentdirs = yaml_stem_init + '/' + shortened.replace('___', '/')
                Path('outdir/' + parentdirs).mkdir(parents=True, exist_ok=True)
                source = 'provenance/workflow/' + location
                # NOTE: Even though we are using subdirectories (not just a single output directory),
                # there is still the possibility of filename collisions, i.e. when scattering.
                # For now, let's use a similar trick as cwltool of append _2, _3 etc.
                # except do it BEFORE the extension.
                # This could still cause problems with slicing, i.e. if you scatter across
                # indices 11-20 first, then 1-10 second, the output file indices will get switched.
                dest = 'outdir/' + parentdirs +  '/' + basename
                if dest in dests:
                    idx = 2
                    while Path(dest).exists():
                        stem = Path(basename).stem
                        suffix = Path(basename).suffix
                        dest = 'outdir/' + parentdirs +  '/' + stem + f'_{idx}' + suffix
                        idx += 1
                dests.add(dest)
                cmd = ['cp', source, dest]
                sub.run(cmd, check=True)


def stage_input_files(yml_inputs: Yaml, root_yml_dir_abs: Path,
                      relative_run_path: bool = True, throw: bool = True) -> None:
    """Copies the input files in yml_inputs to the working directory.

    Args:
        yml_inputs (Yaml): The yml inputs file for the root workflow.
        root_yml_dir_abs (Path): The absolute path of the root workflow yml file.
        relative_run_path (bool): Controls whether to use subdirectories or\n
        just one directory when writing the compiled CWL files to disk
        throw (bool): Controls whether to raise/throw a FileNotFoundError.

    Raises:
        FileNotFoundError: If throw and it any of the input files do not exist.
    """
    for key, val in yml_inputs.items():
        if isinstance(val, Dict) and val.get('class', '') == 'File':
            path = root_yml_dir_abs / Path(val['path'])
            if not path.exists() and throw:
                #raise FileNotFoundError(f'Error! {path} does not exist!')
                print(f'Error! {path} does not exist! (Did you forget to use an explicit edge?)')
                print('See https://workflow-inference-compiler.readthedocs.io/en/latest/userguide.html#explicit-edges')
                sys.exit(1)

            relpath = Path('autogenerated/') if relative_run_path else Path('.')
            pathauto = relpath / Path(val['path']) # .name # NOTE: Use .name ?
            pathauto.parent.mkdir(parents=True, exist_ok=True)

            if path != pathauto:
                cmd = ['cp', str(path), str(pathauto)]
                proc = sub.run(cmd, check=False)


if __name__ == '__main__':
    main()
