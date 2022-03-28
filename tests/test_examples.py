from pathlib import Path
import pytest
import subprocess as sub
import sys
from unittest.mock import patch

import graphviz
import yaml

from wic import auto_gen_header
import wic.cli
import wic.compiler
import wic.main
import wic.utils


@pytest.mark.slow
def test_examples() -> None:
    testargs = ['wic', '--yaml', '', '--cwl_output_intermediate_files', 'True']  # ignore --yaml
    # For now, we need to enable --cwl_output_intermediate_files. See comment in compiler.py
    with patch.object(sys, 'argv', testargs):
        args = wic.cli.parser.parse_args()

    tools_cwl = wic.main.get_tools_cwl(Path('.'))
    yml_paths = wic.main.get_yml_paths(Path('examples/gromacs/'))

    # First compile all of the workflows.
    for yml_path in yml_paths:
        rootgraph = graphviz.Digraph(name=f'cluster_{yml_path}')
        rootgraph.attr(newrank='True')
        workflow_data = wic.compiler.compile_workflow(args, [], [rootgraph], {}, tools_cwl, True, yml_paths[yml_path], yml_paths)
        recursive_data = workflow_data[0]
        sub_node_data = recursive_data[0]
        yaml_stem = sub_node_data[0]

        wic.utils.write_to_disk(recursive_data)

        # Now blindly run all workflows and (if all inputs are present) check for return code 0.
        # Workflows are first validated before runtime, so this also checks for validity.
        # NOTE: Do not use --cachedir; we want to actually test everything.
        cmd = ['cwltool', '--outdir', f'outdir/{yaml_stem}', f'{yaml_stem}.cwl', f'{yaml_stem}_inputs.yml']
        proc = sub.run(cmd, stdout=sub.PIPE, stderr=sub.STDOUT)  # Capture the output
        if not proc.returncode == 0:
            # Since some of the workflows will be subworkflows
            # (i.e. will not have all inputs), we need to check for
            # "Missing required input parameter" and only fail the
            # workflows which should have succeeded.
            missing_input = "Missing required input parameter"
            bad_format = "File has an incompatible format"
            output = proc.stdout.decode("utf-8")
            if not (missing_input in output or bad_format in output):
                print(f"Error! {yml_paths[yml_path]} failed!")
                print(output)
                assert proc.returncode == 0


def test_cwl_embedding_independence() -> None:
    testargs = ['wic', '--yaml', '']  # ignore --yaml
    with patch.object(sys, 'argv', testargs):
        args = wic.cli.parser.parse_args()

    tools_cwl = wic.main.get_tools_cwl(Path('.'))
    yml_paths = wic.main.get_yml_paths(Path('examples/gromacs'))

    # First compile all of the workflows once.
    for yml_path in yml_paths:
        rootgraph = graphviz.Digraph(name=f'cluster_{yml_path}')
        rootgraph.attr(newrank='True')
        workflow_data = wic.compiler.compile_workflow(args, [], [rootgraph], {}, tools_cwl, True, yml_paths[yml_path], yml_paths)
        recursive_data = workflow_data[0]
        node_data_lst = wic.utils.flatten_recursive_data(recursive_data)

        # Now, for each subworkflow of the given root workflow, compile the
        # subworkflow again from scratch, as if it were the root workflow,
        # and check that the generated CWL is identical. In other words,
        # check that the generated CWL of a subworkflow is independent of its
        # embedding into a parent workflow.
        for sub_node_data in node_data_lst:
            (sub_yml_path_str, sub_cwl_embedded, sub_yaml_inputs_, sub_rootgraph_) = sub_node_data

            rootgraph_fakeroot = graphviz.Digraph(name=f'cluster_{sub_yml_path_str}')
            rootgraph_fakeroot.attr(newrank='True')
            workflow_data_fakeroot = wic.compiler.compile_workflow(args, [], [rootgraph_fakeroot], {}, tools_cwl, True, yml_paths[sub_yml_path_str], yml_paths)
            sub_node_data_fakeroot = workflow_data_fakeroot[0][0]
            sub_cwl_fakeroot = sub_node_data_fakeroot[1]

            if not sub_cwl_embedded == sub_cwl_fakeroot:
                yaml_content = yaml.dump(sub_cwl_embedded, sort_keys=False, line_break='\n', indent=2)
                with open(f'{sub_yml_path_str}_embedded.cwl', 'w') as w:
                    w.write('#!/usr/bin/env cwl-runner\n')
                    w.write(auto_gen_header)
                    w.write(''.join(yaml_content))
                yaml_content = yaml.dump(sub_cwl_fakeroot, sort_keys=False, line_break='\n', indent=2)
                with open(f'{sub_yml_path_str}_fakeroot.cwl', 'w') as w:
                    w.write('#!/usr/bin/env cwl-runner\n')
                    w.write(auto_gen_header)
                    w.write(''.join(yaml_content))
            assert sub_cwl_embedded == sub_cwl_fakeroot
