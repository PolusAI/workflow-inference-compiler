import argparse
from pathlib import Path

import subprocess as sub
import sys
from typing import Any, Dict
import zipfile

import yaml


def check_options(options: Dict[str, Any], args: argparse.Namespace) -> None:
    """Performs sanity checks for various settings in the `options` tag.

    Args:
        options (Dict[str, Any]): The `options` tag in the root of the yaml config file.
        args (argparse.Namespace): The command line arguments.
    """
    if 'verbose' not in options:
        options['verbose'] = True #'yes'
    if 'resume_setup' not in options:
        options['resume_setup'] = True #'yes'
    if 'resume_simulation' not in options:
        options['resume_simulation'] = True #'yes'

    if options.get('start_from_trailblaze_samples', '') == 'yes':
        print('Warning: start_from_trailblaze_samples should always be set to no; setting to no')
    options['start_from_trailblaze_samples'] = False #'no' # yes generates NaNs!

    if 'checkpoint_interval' in options and 'switch_phase_interval' in options:
        if options['switch_phase_interval'] < 0 or options['checkpoint_interval'] < 0:
            print('Error! switch_phase_interval and checkpoint_interval should both be >= 0')
            sys.exit(1)
        if options['switch_phase_interval'] < options['checkpoint_interval']:
            print('Error! switch_phase_interval is less than checkpoint_interval!')
            print('This will cause yank to enter an infinite loop!')
            sys.exit(1)
        ratio: float = options['switch_phase_interval'] / options['checkpoint_interval']
        if not ratio.is_integer():
            print('Error! switch_phase_interval should be an integer multiple of checkpoint_interval!')
            sys.exit(1)
    else:
        default_checkpoint_interval = 10
        default_switch_phase_interval = 10
        options['checkpoint_interval'] = default_checkpoint_interval
        options['switch_phase_interval'] = default_switch_phase_interval

    if 'number_of_equilibration_iterations' not in options:
        options['number_of_equilibration_iterations'] = 25
    if args.phase in ['setup_only', 'trailblaze']:
        options['number_of_equilibration_iterations'] = 0

    if 'default_number_of_iterations' not in options:
        options['default_number_of_iterations'] = options['checkpoint_interval']
    if args.phase in ['setup_only', 'trailblaze', 'equilibration']:
        options['default_number_of_iterations'] = 0

def make_molecules_tag(args: argparse.Namespace) -> Dict[str, Any]:
    """Returns the basic tags which describe the molecules in a system.

    Args:
        args (argparse.Namespace): The command line arguments.

    Returns:
        Dict[str, Any]: The basic tags which describe the molecules in a system.
    """
    mol = {
        'receptor_name': {'filepath': args.input_receptor_path},
        'ligand_name': {'filepath': args.input_ligand_path,
                        'antechamber': {'charge_method': 'bcc'}},
    }
    return mol


def make_systems_tag_group1() -> Dict[str, Any]:
    """Returns the basic tags which describe a molecular system.

    Returns:
        Dict[str, Any]: The basic tags which describe a molecular system.
    """
    systems = {'systems': {'system_name': {
        'receptor': 'receptor_name',
        'ligand': 'ligand_name',
        'solvent': 'spce_50mM_acetate_55',
        'leap': {'parameters': ['oldff/leaprc.ff99SBildn', 'leaprc.gaff']},
    }}}
    return systems


def make_systems_tag_group2(complex_top_filename: str, ligand_top_filename: str,
                            args: argparse.Namespace) -> Dict[str, Any]:
    """Returns the basic tags which describe a molecular system.

    Args:
        complex_top_filename (str): The topology filepath for the complex phase.
        ligand_top_filename (str): The topology filepath for the complex phase.
        args (argparse.Namespace): The command line arguments.

    Returns:
        Dict[str, Any]: The basic tags which describe a molecular system.
    """
    gromacs_include_dir = '/miniconda/share/gromacs/top/'

    systems = {'systems': {'system_name': {
        'phase1_path': [complex_top_filename, args.input_complex_crd_path],
        'phase2_path': [ligand_top_filename, args.input_ligand_crd_path],
        'ligand_dsl': 'resname MOL',
        #'solvent_dsl': 'resname WAT',  # optional
        'solvent': 'spce_50mM_acetate_55',
        'gromacs_include_dir': gromacs_include_dir,
    }}}
    return systems


def make_experiment_tag(system_name: str = 'system_name') -> Dict[str, Any]:
    """Returns the basic tags which describe an experiment.

    Args:
        system_name (str, optional): The user-defined name of the system. Defaults to 'system_name'.

    Returns:
        Dict[str, Any]: The basic tags which describe an experiment.
    """
    # TODO: Automatically determine residue indices
    # residue_indices_str = ' or '.join([f'resi {i}' for i in residue_indices])
    exp = {
        'system': system_name,
        'sampler': 'repex',
        'options': {
            'temperature': '302.15*kelvin'},
        'protocol': 'binding-auto',
    }
        # Cannot use restraints on 32 bit systems :(
        #restraint:
        #  #type: FlatBottom
        #  restrained_ligand_atoms: (resname MOL) and (mass > 1.5)
        #  restrained_receptor_atoms: f'({residue_indices_str}) and (mass > 1.5)'
    return exp


def unzip_topology_files(zip_path: str) -> str:
    """Unzips the topology files and returns the filename of the main topology file.

    Args:
        zip_path (str): The path to the topology zip file.

    Returns:
        str: The name of the main topology file.
    """
    filenames = []
    with zipfile.ZipFile(zip_path, 'r') as zip_file:
        filenames = zip_file.namelist()
        # This extracts to a temporary directory which is not writable when using Docker.
        #zip_file.extractall(Path(outdir).parent)
    # Instead, use unzip which extracts in-place (in the same directory)
    cmd = ['unzip', '-n', zip_path] # Use -n to avoid overwriting files.
    sub.run(cmd, check=True)
    top_filenames = [fn for fn in filenames if Path(fn).suffix == '.top']
    if len(top_filenames) != 1:
        print(f'Error! There should be exactly 1 *.top file: {top_filenames}')
        sys.exit(1)
    top_filename = top_filenames[0]
    return top_filename


def cli() -> argparse.ArgumentParser:
    """Returns a parser for the command line arguments.

    Returns:
        argparse.ArgumentParser: A parser for the command line arguments.
    """
    parser = argparse.ArgumentParser()
    phase_choices = ['setup_only', 'trailblaze', 'equilibration', 'production', 'analyze', 'report']
    parser.add_argument('--phase', choices=phase_choices)
    parser.add_argument('--yaml')
    # group 1
    parser.add_argument('--input_receptor_path')
    parser.add_argument('--input_ligand_path')
    # group 2
    parser.add_argument('--input_complex_top_zip_path')
    parser.add_argument('--input_complex_crd_path')
    parser.add_argument('--input_ligand_top_zip_path')
    parser.add_argument('--input_ligand_crd_path')
    # TODO: Check that group2 and setup_only are mutually exclusive
    return parser


def main() -> None:
    """Fills in a yank config yaml template and runs yank."""
    print('yank_wrapper sys.argv', sys.argv)

    parser = cli()
    args = parser.parse_args()

    group1 = bool(args.input_receptor_path or args.input_ligand_path)
    group2 = bool(args.input_complex_top_zip_path or args.input_complex_crd_path or
                  args.input_ligand_top_zip_path or args.input_ligand_crd_path)
    mutex = (group1 and not group2) or (group2 and not group1)
    if not mutex:
        print('Error! The following inputs are mutually exclusive:')
        print('--input_receptor_path --input_ligand_path')
        print('--input_complex_top_zip_path --input_complex_crd_path --input_ligand_top_zip_path --input_ligand_crd_path')
        sys.exit(1)

    # Unzip topology archives
    complex_top_filename = unzip_topology_files(args.input_complex_top_zip_path)
    ligand_top_filename = unzip_topology_files(args.input_ligand_top_zip_path)

    # Create modified yaml file (update paths, number of iterations, perform some checks, etc)
    yaml_path = Path(args.yaml)
    with open(yaml_path, mode='r', encoding='utf-8') as f:
        yaml_dict: Dict[str, Any] = yaml.safe_load(f.read())

    if group1:
        yaml_dict['molecules'] = make_molecules_tag(args)

    if group1:
        systems_tag = make_systems_tag_group1()
    else: # group2
        systems_tag = make_systems_tag_group2(complex_top_filename, ligand_top_filename, args)

    #yaml_dict['systems'] = {**yaml_dict['systems'], **systems_tag['systems']} # merge?
    yaml_dict['systems'] = systems_tag['systems'] # overwrite

    yaml_dict['experiment_name'] = make_experiment_tag()
    #yaml_dict['experiments'] = yaml_dict.get('experiments', []).append('experiment_name') # append?
    yaml_dict['experiments'] = ['experiment_name'] # overwrite

    options: Dict[str, Any] = yaml_dict['options']
    check_options(options, args)

    yaml_filename = f'{yaml_path.name.split(".")[0]}_{args.phase}.yaml'
    with open(yaml_filename, mode='w', encoding='utf-8') as f:
        f.write(yaml.dump(yaml_dict))

    # Run yank
    subcommand = ['script', '-y', yaml_filename] # trailblaze, equilibration, production
    if args.phase == 'setup_only':
        subcommand = ['script', '--setup-only', '-y', yaml_filename]
    if args.phase == 'analyze':
        subcommand = ['analyze', '--store=output/experiments/', '--fulltraj']
    if args.phase == 'report':
        subcommand = ['analyze', 'report', '--store=output/experiments/', '--output=report.ipynb', '--fulltraj']

    cmd = ['yank'] + subcommand
    print('Running', cmd)
    sub.run(cmd, check=True)


if __name__ == '__main__':
    main()
