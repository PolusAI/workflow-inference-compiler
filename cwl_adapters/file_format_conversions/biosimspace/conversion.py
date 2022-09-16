import argparse
import os
import sys

import BioSimSpace as BSS

# Simplified version of https://github.com/michellab/BioSimSpace/blob/devel/demo/conversion.py
# Supported formats are currently: ['Gro87', 'GroTop', 'MOL2', 'PDB', 'PRM7', 'PSF', 'RST', 'RST7']
# Note that the individual CWL scripts are specialized to specific conversions so that
# we can use more precise edam formats and thus let CWL catch more errors.
# In other words, the cwl scripts are type-safe wrappers around this script.

def main() -> None:
    parser = argparse.ArgumentParser()
    required_args = parser.add_argument_group('required arguments')
    required_args.add_argument('--input_top_path', required=True)
    required_args.add_argument('--input_crd_path', required=True)
    required_args.add_argument('--output_top_format', required=True)
    required_args.add_argument('--output_crd_format', required=True)
    required_args.add_argument('--filename_prefix', required=False, default='top')
    args = parser.parse_args()

    if args.output_top_format not in BSS.IO.fileFormats():
        print(f'Error! Unsupported topology format {args.output_top_format}')
        print("Supported formats are ", BSS.IO.fileFormats())
        sys.exit(1)
    if args.output_crd_format not in BSS.IO.fileFormats():
        print(f'Error! Unsupported coordinate format {args.output_crd_format}')
        print("Supported formats are ", BSS.IO.fileFormats())
        sys.exit(1)

    input_filenames = [args.input_top_path, args.input_crd_path]
    # This path is w.r.t. the Docker image jakefennick/biosimspace
    properties = {'GROMACS_PATH': '/miniconda/share/gromacs/top/'}
    system = BSS.IO.readMolecules(input_filenames, properties)

    output_file_formats = [args.output_top_format, args.output_crd_format]
    prefix = args.filename_prefix
    output_filenames = BSS.IO.saveMolecules(prefix, system, output_file_formats)
    print(output_filenames)

    # Rename output files to use standard extensions. Support more than one filename prefix?
    if os.path.exists(f'{prefix}.prm7'):
        os.rename(f'{prefix}.prm7', f'{prefix}.prmtop')
    if os.path.exists(f'{prefix}.rst7'):
        os.rename(f'{prefix}.rst7', f'{prefix}.inpcrd')
    # BioSimSpace automatically renames *.GroTop to *.top
    # BioSimSpace automatically renames *.Gro87 to *.gro


if __name__ == '__main__':
    main()
