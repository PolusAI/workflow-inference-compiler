import argparse

parser = argparse.ArgumentParser(prog='main', description='Convert a high-level yaml workflow file to CWL.')
parser.add_argument('--yaml', type=str, required=True,
                    help='Yaml workflow file')

parser.add_argument('--cwl_dir', type=str, required=False, default='biobb',
                    help='Directory which contains the CWL CommandLineTools and/or Workflows')
# Change default to True for now. See comment in compiler.py
parser.add_argument('--cwl_output_intermediate_files', type=bool, required=False, default=True,
                    help='Enable output files which are used between steps (for debugging).')
parser.add_argument('--cwl_run_local', type=bool, required=False, default=False,
                    help='After generating the cwl file(s), run it on localhost.')
# NOTE: If cwl_run_slurm is enabled, you MUST enable cwl_inline_subworkflows!
# Plugins with 'class: Workflow' (i.e. subworkflows) are not currently supported.
parser.add_argument('--cwl_run_slurm', type=bool, required=False, default=False,
                    help='After generating the cwl file, run it on labshare using the slurm driver.')
parser.add_argument('--cwl_inline_subworkflows', type=bool, required=False, default=False,
                    help='Before generating the cwl file, inline all subworkflows.')
parser.add_argument('--cwl_validate', type=bool, required=False, default=False,
                    help='After generating the cwl file, validate it.')

parser.add_argument('--compute_url', type=str, required=False, default=False,
                    help='The URL associated with the labshare slurm driver.')

parser.add_argument('--graph_label_edges', type=bool, required=False, default=False,
                    help='Label the graph edges with the name of the intermediate input/output.')
parser.add_argument('--graph_label_stepname', type=bool, required=False, default=False,
                    help='Prepend the step name to each step node.')
parser.add_argument('--graph_show_inputs', type=bool, required=False, default=False,
                    help='Add nodes to the graph representing the workflow inputs.')
parser.add_argument('--graph_show_outputs', type=bool, required=False, default=False,
                    help='Add nodes to the graph representing the workflow outputs.')
parser.add_argument('--graph_inline_depth', type=int, required=False, default=100,  # 100 == all
                    help='Controls the depth of subgraphs which are displayed separately or positioned within the main graph.')