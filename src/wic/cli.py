import argparse
import sys

parser = argparse.ArgumentParser(prog='main', description='Convert a high-level yaml workflow file to CWL.')
parser.add_argument('--yaml', type=str, required=('--generate_schemas_only' not in sys.argv),
                    help='Yaml workflow file')

parser.add_argument('--generate_schemas_only', type=bool, required=False,
                    help='Generate schemas for the files in --cwl_dirs_file and --yml_dirs_file.')
parser.add_argument('--cwl_dirs_file', type=str, required=False, default='cwl_dirs.txt',
                    help='Configuration file which lists the directories which contains the CWL CommandLineTools')
parser.add_argument('--yml_dirs_file', type=str, required=False, default='yml_dirs.txt',
                    help='Configuration file which lists the directories which contains the YAML Workflows')
# Change default to True for now. See comment in compiler.py
parser.add_argument('--cwl_output_intermediate_files', type=bool, required=False, default=True,
                    help='Enable output files which are used between steps (for debugging).')

parser.add_argument('--parallel', type=bool, required=False, default=False,
                    help='''When running locally, execute independent steps in parallel.
                    \nThis is required for real-time analysis, but it may cause issues with
                    \nhanging (particularly when scattering). See user guide for details.''')

group_run = parser.add_mutually_exclusive_group()
group_run.add_argument('--cwl_run_local', type=bool, required=False, default=False,
                    help='After generating the cwl file(s), run it on localhost.')
group_run.add_argument('--cwl_run_slurm', type=bool, required=False, default=False,
                    help='After generating the cwl file, run it on labshare using the slurm driver.')
# Use required=('--cwl_run_slurm' in sys.argv) make other args conditionally required.
# See https://stackoverflow.com/questions/19414060/argparse-required-argument-y-if-x-is-present
# For example, if cwl_run_slurm is enabled, you MUST enable cwl_inline_subworkflows!
# Plugins with 'class: Workflow' (i.e. subworkflows) are not currently supported.

parser.add_argument('--cwl_inline_subworkflows', type=bool, required=False, default=('--cwl_run_slurm' in sys.argv),
                    help='Before generating the cwl file, inline all subworkflows. Required for --cwl_run_slurm')
parser.add_argument('--cwl_inference_use_naming_conventions', type=bool, required=False, default=False,
                    help='Enables the use of naming conventions in the inference algorithm')
parser.add_argument('--cwl_validate', type=bool, required=False, default=False,
                    help='After generating the cwl file, validate it.')
parser.add_argument('--cachedir', type=str, required=False, default='cachedir',
                    help='The directory to save intermediate results; useful with RealtimePlots.py')

aws_url = 'http://compute.ci.aws.labshare.org'
ncats_url = 'https://compute.scb-ncats.io/'

parser.add_argument('--compute_url', type=str, default=ncats_url,
                    help='The URL associated with the labshare slurm driver. Required for --cwl_run_slurm')
parser.add_argument('--compute_access_token', type=str, required=('--cwl_run_slurm' in sys.argv),
                    help="""The access_token used for authentication. Required for --cwl_run_slurm
                    For now, get this manually from https://a-qa.labshare.org/""")

parser.add_argument('--graph_label_edges', type=bool, required=False, default=False,
                    help='Label the graph edges with the name of the intermediate input/output.')
parser.add_argument('--graph_label_stepname', type=bool, required=False, default=False,
                    help='Prepend the step name to each step node.')
parser.add_argument('--graph_show_inputs', type=bool, required=False, default=False,
                    help='Add nodes to the graph representing the workflow inputs.')
parser.add_argument('--graph_show_outputs', type=bool, required=False, default=False,
                    help='Add nodes to the graph representing the workflow outputs.')
parser.add_argument('--graph_inline_depth', type=int, required=False, default=sys.maxsize,
                    help='Controls the depth of subgraphs which are displayed.')
parser.add_argument('--graph_dark_theme', type=bool, default=False,
                    help='Changees the color of the fonts and edges from white to black.')
