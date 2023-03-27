from typing import Dict, List
from ..wic_types import Json
from .gromacs_mdp import gromacs_mdp_schema, gromacs_selection_groups

# The following schema are used to provide IntelliSense code completion for
# the biob config: tags. The CWL tools that appear in the examples/ have been
# manually curated below (to add detailed enums, etc). Otherwise, we could also
# attempt to screen-scrape the associated URLs as in gromacs_mdp.py


def default_schema(url: bool = False) -> Json:
    """A basic default schema (to avoid copy & paste).

    Args:
        url (bool, optional): Determines whether to include the $schema url. Defaults to False.

    Returns:
        Json: A basic default schema
    """
    schema: Json = {}
    schema['type'] = 'object'
    schema['additionalProperties'] = False
    if url:
        schema['$schema'] = 'https://json-schema.org/draft/2020-12/schema'
    return schema


def biobb_container_schema() -> Json:
    # pylint: disable=line-too-long
    """A schema for common container-related options.

    Returns:
        Json: A schema for common container-related options.
    """
    gmx_lib = {'type': 'string', 'description': '(None) Path set GROMACS GMXLIB environment variable.'}
    gmx_path = {'type': 'string', 'description': '(“gmx”) Path to the GROMACS executable binary.'}
    remove_tmp = {'type': 'boolean', 'description': '(True) [WF property] Remove temporal files.'}
    restart = {'type': 'boolean', 'description': '(False) [WF property] Do not execute if output files exist.'}
    container_path = {'type': 'string', 'description': '(None) Path to the binary executable of your container.'}
    container_image = {'type': 'string', 'description': '(“gromacs/gromacs:latest”) Container Image identifier.'}
    container_volume_path = {'type': 'string',
                             'description': '(“/data”) Path to an internal directory in the container.'}
    container_working_dir = {'type': 'string', 'description': '(None) Path to the internal CWD in the container.'}
    container_user_id = {'type': 'string', 'description': '(None) User number id to be mapped inside the container.'}
    container_shell_path = {'type': 'string',
                            'description': '(“/bin/bash”) Path to the binary executable of the container shell.'}
    container = {'gmx_lib': gmx_lib, 'gmx_path': gmx_path,
                 'remove_tmp': remove_tmp, 'restart': restart,
                 'container_path': container_path,
                 'container_image': container_image,
                 'container_volume_path': container_volume_path,
                 'container_working_dir': container_working_dir,
                 'container_user_id': container_user_id,
                 'container_shell_path': container_shell_path}
    return container


def biobb_selection_schema() -> Json:
    """A schema for gromacs selection groups.

    Returns:
        Json: A schema for gromacs selection groups.
    """
    # See https://biobb-analysis.readthedocs.io/en/latest/gromacs.html#module-gromacs.gmx_rms
    # See https://biobb-analysis.readthedocs.io/en/latest/gromacs.html#module-gromacs.gmx_rgyr
    xvg = {'type': 'string', 'enum': ['xmgrace', 'xmgr', 'none'], 'description': '(“none”) XVG plot formatting.'}

    groups_schemas = []
    for group, desc in gromacs_selection_groups().items():
        g_s = {'type': 'string', 'const': group, 'description': desc}
        groups_schemas.append(g_s)
    groups_schemas.append({'type': 'string', 'pattern': '^resname.*'})  # TODO: Check for other patterns
    # NOTE: Use anyOf instead of enum so we can add descriptions to each value
    desc = """(“System”) Group where the rms will be performed.
    If input_index_path provided, check the file for the accepted values."""
    groups = {'anyOf': groups_schemas, 'description': desc}
    # groups = {'type': 'string', 'enum': [s['const'] for s in groups_schemas]}
    # NOTE: This schema appears to trigger some kind of a bug in the vscode YAML
    # extension. Specifically, IntelliSense code completion will not work when
    # there is not enough information for at least one exact match. For example,
    # if you type "MainChai" and press ctrl-space, nothing happens. But once you
    # #type "MainChain", IntelliSense correctly suggests "MainChain", "MainChain+Cb",
    # and "MainChain+H".

    schema = default_schema()
    schema['properties'] = {'xvg': xvg, 'selection': groups,
                            **biobb_container_schema()}  # no gmx_lib ??
    return schema


def biobb_gmx_energy_schema() -> Json:
    """A schema for gromacs gmx_energy options.

    Returns:
        Json: A schema for gromacs gmx_energy options.
    """
    # See https://biobb-analysis.readthedocs.io/en/latest/gromacs.html#module-gromacs.gmx_energy
    xvg = {'type': 'string', 'enum': ['xmgrace', 'xmgr', 'none'], 'description': '(“none”) XVG plot formatting.'}
    # TODO: Check this list; the documentation appears to have typos...
    terms_list = ['Angle', 'Proper-Dih.', 'Improper-Dih.', 'LJ-14', 'Coulomb-14',
                  'LJ-(SR)', 'Coulomb-(SR)', 'Coul.-recip.', 'Position-Rest.', 'Potential',
                  'Kinetic-En.', 'Total-Energy', 'Temperature', 'Pressure', 'Constr.-rmsd',
                  'Box-X', 'Box-Y', 'Box-Z', 'Volume', 'Density', 'pV', 'Enthalpy', 'Vir-XX',
                  'Vir-XY', 'Vir-XZ', 'Vir-YX', 'Vir-YY', 'Vir-YZ', 'Vir-ZX', 'Vir-ZY',
                  'Vir-ZZ', 'Pres-XX', 'Pres-XY', 'Pres-XZ', 'Pres-YX', 'Pres-YY',
                  'Pres-YZ', 'Pres-ZX', 'Pres-ZY', 'Pres-ZZ', 'Surf*SurfTen', 'Box-Vel-XX',
                  'Box-Vel-YY', 'Box-Vel-ZZ', 'Mu-X', 'Mu-Y', 'Mu-Z', 'T-Protein',
                  'T-non-Protein', 'Lamb-Protein', 'Lamb-non-Protein']
    terms_schema = {'type': 'string', 'enum': terms_list}
    terms = {'type': 'array', 'items': terms_schema, 'description': '([“Potential”]) Energy terms.'}

    schema = default_schema()
    schema['properties'] = {'xvg': xvg, 'terms': terms,
                            **biobb_container_schema()}  # no gmx_lib ??
    return schema


def biobb_genion_schema() -> Json:
    """A schema for gromacs genion options.

    Returns:
        Json: A schema for gromacs genion options.
    """
    # See https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.genion
    replaced_group = {'type': 'string', 'description':
                      '(“SOL”) Group of molecules that will be replaced by the solvent.'}
    neutral = {'type': 'boolean', 'description': '(False) Neutralize the charge of the system.'}
    concentration = {'type': 'number', 'description': '(0.05) [0~10|0.01] Concentration of the ions in (mol/liter).'}
    seed = {'type': 'number', 'description': '(1993) Seed for random number generator.'}

    schema = default_schema()
    schema['properties'] = {'replaced_group': replaced_group, 'neutral': neutral,
                            'concentration': concentration, 'seed': seed,
                            **biobb_container_schema()}
    return schema


def biobb_solvate_schema() -> Json:
    """A schema for gromacs solvate options.

    Returns:
        Json: A schema for gromacs solvate options.
    """
    # See https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.solvate
    desc = '(0.0) [0~100|0.1] Thickness in nanometers of optional water layer around solute.'
    shell = {'type': 'number', 'description': desc}

    schema = default_schema()
    schema['properties'] = {'shell': shell,
                            **biobb_container_schema()}
    return schema


def biobb_editconf_schema() -> Json:
    """A schema for gromacs editconf options.

    Returns:
        Json: A schema for gromacs editconf options.
    """
    # See https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.editconf
    desc = '(1.0) [0~100|0.1] Distance of the box from the outermost atom in nm. ie 1.0nm = 10 Angstroms.'
    distance_to_molecule = {'type': 'number', 'description': desc}
    box_vals = ['cubic', 'triclinic', 'dodecahedron', 'octahedron']
    desc = """(“cubic”) Geometrical shape of the solvent box. Values:
    cubic (rectangular box with all sides equal),
    triclinic (triclinic box),
    dodecahedron (rhombic dodecahedron),
    octahedron (truncated octahedron)."""
    box_type = {'type': 'string', 'enum': box_vals, 'description': desc}
    center_molecule = {'type': 'boolean', 'description': '(True) Center molecule in the box.'}
    box_vector_lengths = {'type': 'array', 'items': {'type': 'number'},
                          'description': '(None) Array of floats defining the box vector lengths ie "0.5 0.5 0.5".'
                          'If this option is used the distance_to_molecule property will be ignored.'}

    schema = default_schema()
    schema['properties'] = {'distance_to_molecule': distance_to_molecule,
                            'box_type': box_type, 'center_molecule': center_molecule, 'box_vector_lengths': box_vector_lengths,
                            **biobb_container_schema()}
    return schema


def biobb_mdrun_schema() -> Json:
    # pylint: disable=line-too-long
    """A schema for gromacs mdrun options.

    Returns:
        Json: A schema for gromacs mdrun options.
    """
    # See https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.mdrun
    mpi_bin = {'type': 'string', 'description': '(None) Path to the MPI runner. Usually “mpirun” or “srun”.'}
    mpi_np = {'type': 'number',
              'description': '(0) [0~1000|1] Number of MPI processes. Usually an integer bigger than 1.'}
    mpi_flags = {'type': 'string', 'description': '(None) Path to the MPI hostlist file.'}
    checkpoint_time = {
        'type': 'number', 'description': '(15) [0~1000|1] Checkpoint writing interval in minutes. Only enabled if an output_cpt_path is provided.'}
    num_threads = {'type': 'number',
                   'description': '(0) [0~1000|1] Let GROMACS guess. The number of threads that are going to be used.'}
    num_threads_mpi = {
        'type': 'number', 'description': '(0) [0~1000|1] Let GROMACS guess. The number of GROMACS MPI threads that are going to be used.'}
    num_threads_omp = {
        'type': 'number', 'description': '(0) [0~1000|1] Let GROMACS guess. The number of GROMACS OPENMP threads that are going to be used.'}
    num_threads_omp_pme = {
        'type': 'number', 'description': '(0) [0~1000|1] Let GROMACS guess. The number of GROMACS OPENMP_PME threads that are going to be used.'}
    use_gpu = {'type': 'boolean', 'description': '(False) Use settings appropriate for GPU. Adds: -nb gpu -pme gpu'}
    gpu_id = {'type': 'string', 'description': '(None) List of unique GPU device IDs available to use.'}
    gpu_tasks = {'type': 'string',
                 'description': '(None) List of GPU device IDs, mapping each PP task on each node to a device.'}

    schema = default_schema()
    schema['properties'] = {'mpi_bin': mpi_bin, 'mpi_np': mpi_np,
                            'mpi_flags': mpi_flags, 'checkpoint_time': checkpoint_time,
                            'num_threads': num_threads, 'num_threads_mpi': num_threads_mpi,
                            'num_threads_omp': num_threads_omp, 'num_threads_omp_pme': num_threads_omp_pme,
                            'use_gpu': use_gpu, 'gpu_id': gpu_id, 'gpu_tasks': gpu_tasks,
                            **biobb_container_schema()}
    return schema


def biobb_pdb_schema() -> Json:
    """A schema for downloading pdb files.

    Returns:
        Json: A schema for downloading pdb files.
    """
    # See https://biobb-io.readthedocs.io/en/latest/api.html#module-api.pdb
    pdb_code = {'type': 'string', 'minLength': 4, 'maxLength': 4, 'description': 'RCSB PDB code'}
    api_id = {'type': 'string', 'enum': ['pdbe', 'pdb', 'mmb'], 'description': 'PDB REST API mirror'}
    schema = default_schema()
    schema['properties'] = {'pdb_code': pdb_code, 'api_id': api_id}
    return schema


def biobb_pdb2gmx_schema() -> Json:
    """A schema for gromacs pdb2gmx options.

    Returns:
        Json: A schema for gromacs pdb2gmx options.
    """
    # See https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.pdb2gmx
    watervals = ['spc', 'spce', 'tip3p', 'tip4p', 'tip5p', 'tips3p']
    water_type = {'type': 'string', 'enum': watervals, 'description': '(“spce”) Water molecule type.'}
    ffvals = ['gromos45a3', 'charmm27', 'gromos53a6', 'amber96', 'amber99',
              'gromos43a2', 'gromos54a7', 'gromos43a1', 'amberGS', 'gromos53a5',
              'amber99sb', 'amber03', 'amber99sb-ildn', 'oplsaa', 'amber94',
              'amber99sb-star-ildn-mut']
    desc = '(“amber99sb-ildn”) Force field to be used during the conversion.'
    forcefield = {'type': 'string', 'enum': ffvals, 'description': desc}
    ignh = {'type': 'boolean', 'description': '(False) Should pdb2gmx ignore the hydrogens in the original structure.'}
    his = {'type': 'string', 'description': '(None) Histidine protonation array.'}
    merge = {'type': 'boolean', 'description': '(False) Merge all chains into a single molecule.'}

    schema = default_schema()
    schema['properties'] = {'water_type': water_type, 'force_field': forcefield,
                            'ignh': ignh, 'his': his, 'merge': merge,
                            **biobb_container_schema()}
    return schema


def biobb_grompp_schema() -> Json:
    """A schema for gromacs grompp options.

    Returns:
        Json: A schema for gromacs grompp options.
    """
    # See https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.grompp
    mdp = default_schema()
    mdp['properties'] = gromacs_mdp_schema()

    schema = default_schema(url=True)
    schema['$id'] = 'gromacs_mdp'
    schema['properties'] = {'mdp': mdp, 'maxwarn': {'type': 'number', "minimum": 0}}
    return schema


def biobb_str_check_add_hydrogens_schema() -> Json:
    """A schema for str_check_add_hydrogens options.

    Returns:
        Json: A schema for str_check_add_hydrogens options.
    """
    # See https://biobb-structure-utils.readthedocs.io/en/latest/utils.html#utils-str-check-add-hydrogens-module
    charges = {'type': 'boolean', 'description':
               '(False) Whether or not to add charges to the output file. If True the output is in PDBQT format.'}
    mode = {'type': 'string', 'enum': ['auto', 'list', 'ph'], 'description': '(auto) Selection mode'}
    p_h = {'type': 'number', 'minimum': 0, 'maximum': 14, 'description':
           '(7.4) Add hydrogens appropriate for pH. Only in case mode ph selected.'}
    lst = {'type': 'string', 'description':
           'List of residues to modify separated by commas (i.e HISA234HID,HISB33HIE). Only in case mode list selected.'}
    keep_canonical_resnames = {'type': 'boolean', 'description': '(False) Whether or not keep canonical residue names'}
    binary_path = {'type': 'string', 'description': '(“check_structure”) path to the check_structure application'}
    remove_tmp = {'type': 'boolean', 'description': '(True) Remove temporal files.'}
    restart = {'type': 'boolean', 'description': '(False) Do not execute if output files exist.'}

    schema = default_schema(url=True)
    schema['properties'] = {'charges': charges, 'mode': mode, 'ph': p_h, 'list': lst,
                            'keep_canonical_resnames': keep_canonical_resnames,
                            'binary_path': binary_path, 'remove_tmp': remove_tmp, 'restart': restart}
    return schema


def biobb_pdb4amber_run_schema() -> Json:
    """A schema for pdb4amber_run options.

    Returns:
        Json: A schema for pdb4amber_run options.
    """
    # See https://biobb-amber.readthedocs.io/en/latest/pdb4amber.html#module-pdb4amber.pdb4amber_run
    remove_hydrogens = {'type': 'boolean', 'description':
                        '(False) Remove hydrogen atoms from the PDB file.'}
    remove_waters = {'type': 'boolean', 'description':
                     '(False) Remove water molecules from the PDB file.'}
    constant_ph = {'type': 'boolean', 'description':
                   '(False) Rename ionizable residues e.g. GLU,ASP,HIS for constant pH simulation.'}
    binary_path = {'type': 'string', 'description': '(“pdb4amber”) Path to the pdb4amber executable binary.'}
    remove_tmp = {'type': 'boolean', 'description': '(True) Remove temporal files.'}
    restart = {'type': 'boolean', 'description': '(False) Do not execute if output files exist.'}

    schema = default_schema(url=True)
    schema['properties'] = {'remove_hydrogens': remove_hydrogens, 'remove_waters': remove_waters,
                            'constant_pH': constant_ph,  # NOTE: capital H
                            'binary_path': binary_path, 'remove_tmp': remove_tmp, 'restart': restart}
    return schema


def cpptraj_mask() -> List[Dict[str, str]]:
    """ Mask options for the cpptraj functions

    Returns:
        List[Dict[str, str]]: A list of dicts
    """
    names = ['all-atoms', 'c-alpha', 'backbone', 'heavy-atoms', 'side-chain', 'solute', 'ions', 'solvent']
    descs = ['All system atoms', 'All c-alpha atoms; protein only', 'Backbone atoms', 'All not backbone atoms',
             'All system atoms except solvent atoms', 'All ion molecules', 'All solvent atoms']
    masks = [{'type': 'string', 'const': name, 'description': desc} for name, desc in zip(names, descs)]
    return masks


def cpptraj_reference() -> List[Dict[str, str]]:
    """ Reference options for the cpptraj functions

    Returns:
         List[Dict[str, str]]: A list of dicts
    """

    names = ['firts', 'average', 'experimental']
    descs = ['Use the first trajectory frame as reference',
             'Use the average of all trajectory frames as reference',
             'Use the experimental structure as reference']

    refs = [{'type': 'string', 'const': name, 'description': desc} for name, desc in zip(names, descs)]
    return refs


def cpptraj_fields() -> List[Dict[str, str]]:
    """ Fields options for cpptraj function

    Returns:
         List[Dict[str, str]]: A list of dicts
    """
    range = '[1~100000|1]'
    start = {'type': 'number', f'description': '(1) {range} Starting frame for slicing'}
    end = {'type': 'number', f'description': '(-1) {range} Ending frame for slicing'}
    steps = {'type': 'number', f'description': '(1) {range} Step for slicing'}

    return [start, end, steps]


def cpptraj_fileformat() -> List[Dict[str, str]]:
    """ File format options for the cpptraj functions

    Returns:
        List[Dict[str, str]]: A list of dicts
    """
    names = ['crd', 'cdf', 'netcdf', 'nc', 'restart', 'ncrestart', 'restartnc', 'dcd', 'charm'
             'cor', 'pdb', 'mol2', 'trr', 'gro', 'binpos', 'xtc', 'cif', 'arc', 'sqm', 'sdf', 'conflib']
    netcdf = 'Format used by netCDF software library for writing and reading chromatography-MS data files'
    amber = 'AMBER coordinate/restart file with 6 coordinates per line'
    descs = ['AMBER trajectory format',
             netcdf, netcdf, netcdf, amber, amber, amber,
             'AMBER trajectory format',
             'Format of CHARMM Residue Topology Files (RTF)', 'Charmm COR',
             'Protein Data Bank format',  'Complete and portable representation of a SYBYL molecule',
             'Trajectory of a simulation experiment used by GROMACS',
             'Translation of the ASCII atom coordinate format to binary code',
             'Portable binary format for trajectories produced by GROMACS package',
             'Entry format of PDB database in mmCIF format', 'Tinker ARC',
             'Tinker ARC', 'SQM Input',
             'One of a family of chemical-data file formats developed by MDL Information Systems',
             'LMOD Conflib']

    file_format = [{'type': 'string', 'const': name, 'description': desc} for name, desc in zip(names, descs)]
    return file_format


def biobb_cpptraj_rms_schema() -> Json:
    """A schema for cpptraj_rms options.

    Returns:
        Json: A schema for cpptraj_rms options.
    """
    # See https://biobb-analysis.readthedocs.io/en/latest/ambertools.html#module-ambertools.cpptraj_rms

    start, end, steps = cpptraj_fields()
    mask = {'anyOf': cpptraj_mask(),
            'description': 'Mask definition'}
    reference = {'anyOf': cpptraj_reference(),
                 'description': 'Reference definition'}
    schema = default_schema()
    schema['properties'] = {'start': start, 'end': end, 'steps': steps,
                            'mask': mask, 'reference': reference}
    return schema


def biobb_cpptraj_rms_nofit_schema() -> Json:
    """A schema for cpptraj_rms_nofit options.

    Returns:
        Json: A schema for cpptraj_rms_nofit options.
    """
    # See https://biobb-analysis.readthedocs.io/en/latest/ambertools.html#module-ambertools.cpptraj_rms

    start, end, steps = cpptraj_fields()
    mask = {'anyOf': cpptraj_mask(),
            'description': 'Mask definition'}
    reference = {'anyOf': cpptraj_reference(),
                 'description': 'Reference definition'}
    nofit = {'type': 'boolean', 'description': '(False) Do not perform best-fit RMSD.'}
    norotate = {'type': 'boolean',
                'description': '(False) If calculating best-fit RMSD, translate but do not rotate coordinates.'}
    nomod = {'type': 'boolean', 'description': '(False) If calculating best-fit RMSD, do not modify coordinates.'}
    schema = default_schema()
    schema['properties'] = {'start': start, 'end': end, 'steps': steps,
                            'mask': mask, 'reference': reference,
                            'nofit': nofit, 'norotate': norotate, 'nomod': nomod}
    return schema


def biobb_cpptraj_rgyr_schema() -> Json:
    """A schema for biobb_cpptraj_rgyr options.

    Returns:
        Json: A schema for biobb_cpptraj_rgyr options.
    """
    # See https://biobb-analysis.readthedocs.io/en/latest/ambertools.html#module-ambertools.cpptraj_rgyr
    start, end, steps = cpptraj_fields()
    mask = {'anyOf': cpptraj_mask(),
            'description': 'Mask definition'}
    schema = default_schema()
    schema['properties'] = {'start': start, 'end': end, 'steps': steps, 'mask': mask}
    return schema


def cpptraj_selection_schema() -> Json:
    """A schema for cpptraj_convert options.

    Returns:
        Json: A schema for cpptraj_convert options.
    """
    # See https://biobb-analysis.readthedocs.io/en/latest/ambertools.html#module-ambertools.cpptraj_convert
    # See https://biobb-analysis.readthedocs.io/en/latest/ambertools.html#module-ambertools.cpptraj_image

    start, end, steps = cpptraj_fields()
    mask = {'anyOf': cpptraj_mask(),
            'description': 'Mask definition'}
    file_format = {'anyOf': cpptraj_fileformat(),
                   'description': 'Output trajectory format'}
    schema = default_schema()
    schema['properties'] = {'start': start, 'end': end, 'steps': steps, 'mask': mask, 'format': file_format}
    return schema


def cpptraj_strip_ambmask_schema() -> Json:
    """A schema for cpptraj_strip with Amaber mask options.

    Returns:
        Json: A schema for cpptraj_strip with Amaber mask options.
    """
    # See https://biobb-analysis.readthedocs.io/en/latest/ambertools.html#module-ambertools.cpptraj_strip
    # https://biobb-analysis.readthedocs.io/en/latest/ambertools.html#module-ambertools.cpptraj_image

    start, end, steps = cpptraj_fields()
    mask = {'type': 'string',  'description': 'Amber atom selection mask definition'}
    file_format = {'anyOf': cpptraj_fileformat(),
                   'description': 'Output trajectory format'}
    schema = default_schema()
    schema['properties'] = {'start': start, 'end': end, 'steps': steps, 'mask': mask, 'format': file_format}
    return schema


def make_ndx_schema() -> Json:
    """A schema for make_ndx options.

    Returns:
        Json: A schema for make_ndx options.
    """
    # See https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.make_ndx
    selection = {'type': 'string', 'description': '"a CA C N O" Heavy atoms. Atom selection string'}

    schema = default_schema()
    schema['properties'] = {'selection': selection}
    return schema


def genrestr_schema() -> Json:
    """A schema for genrestr options.

    Returns:
        Json: A schema for genrestr options.
    """
    # See https://biobb-md.readthedocs.io/en/latest/gromacs.html?highlight=genrestr#module-gromacs.genrestr

    restrained_group = {'type': 'string', 'description': 'Index group that will be restrained.'}
    force_constants = {'type': 'string', 'description': 'Array of three floats defining the force constants'}
    schema = default_schema()
    schema['properties'] = {'restrained_group': restrained_group, 'force_constants': force_constants}
    return schema


def gmx_trjconv_str_schema() -> Json:
    """A schema for gmx_trjconv_str options.

    Returns:
        Json: A schema for gmx_trjconv_str options.
    """
    # See https://biobb-analysis.readthedocs.io/en/latest/gromacs.html#module-gromacs.gmx_trjconv_str

    selection = {'type': 'string', 'description': 'Group where the trjconv will be performed'}
    skip = {'type': 'number', 'description': '[0~10000|1] Only write every nr-th frame.'}
    start = {'type': 'number', 'description':
             '[0~10000|1] Time of first frame to read from trajectory (default unit ps).'}
    end = {'type': 'number', 'description':
           '[0~10000|1] Time of last frame to read from trajectory (default unit ps).'}
    dt = {'type': 'number', 'description':
          'Only write frame when t MOD dt = first time (ps).'}
    output_name = {'type': 'string', 'description': ' File name for ensemble of output files'}
    output_type = {'type': 'string', 'description': 'File type for ensemble of output files.'}

    schema = default_schema()
    schema['properties'] = {'selection': selection, 'skip': skip, 'start': start, 'end': end,
                            'dt': dt, 'output_name': output_name, 'output_type': output_type}
    return schema


def leap_gen_top_schema() -> Json:
    """A schema for leap_gen_top options.

    Returns:
        Json: A schema for leap_gen_top options.
    """
    # See https://biobb-amber.readthedocs.io/en/latest/leap.html#module-leap.leap_gen_top

    forcefield = {'type': 'array', 'items': {'type': 'string'},
                  'description': 'Forcefield to be used for the structure generation'}
    schema = default_schema()
    schema['properties'] = {'forcefield': forcefield}
    return schema


def biobb_process_minout_schema() -> Json:
    """A schema for process_minout options.

    Returns:
        Json: A schema for process_minout options.
    """
    # See https://biobb-amber.readthedocs.io/en/latest/process.html#module-process.process_minout

    minout_energy_terms = ['ANGLE', 'BOND', 'DIHEDRAL', 'EEL', 'EEL14',
                           'ENERGY', 'GMAX', 'HBOND', 'NAME', 'NSTEP', 'NUMBER',
                           'RESTRAINT', 'RMS', 'VDW14', 'VDWAALS']

    minout_term_amber = {'type': 'string', 'enum': minout_energy_terms}
    minout_terms_amber = {'type': 'array', 'items': minout_term_amber,
                          'description': '([“ENERGY”]) Statistics descriptors.'}
    return {**default_schema(), 'properties': {'terms': minout_terms_amber}}


def biobb_leap_solvate_schema() -> Json:
    """A schema for leap_solvate options.

    Returns:
        Json: A schema for leap_solvate options.
    """
    # See https://biobb-amber.readthedocs.io/en/latest/leap.html#module-leap.leap_solvate
    forcefield_terms = ['protein.ff14SB', 'protein.ff19SB', 'DNA.bsc1', 'DNA.OL15', 'RNA.OL3', 'gaff']
    forcefield = {'type': 'array', 'items': {'type': 'string', 'enum': forcefield_terms},
                  'description': 'Forcefield to be used for the structure generation.'}
    water_type_terms = ['POL3BOX', 'QSPCFWBOX', 'SPCBOX', 'SPCFWBOX', 'TIP3PBOX', 'TIP3PFBOX',
                        'TIP4PBOX', 'TIP4PEWBOX', 'OPCBOX', 'OPC3BOX', 'TIP5PBOX']
    water_type = {'type': 'string', 'enum': water_type_terms,
                  'description': 'Water molecule parameters to be used for the topology'}
    box_type = {'type': 'string', 'enum': ['cubic', 'truncated_octahedron'],
                'description': '(“truncated_octahedron”) Type for the MD system box'}
    ions_type_terms = ['ionsjc_tip3p', 'ionsjc_spce', 'ionsff99_tip3p', 'ions_charmm22', 'ionsjc_tip4pew', 'None']
    ions_type = {'type': 'string', 'enum': ions_type_terms,
                 'description': '(“ionsjc_tip3p”) Ions type.'}
    neutralise = {'type': 'boolean', 'description':
                  '(“False”) Energetically neutralise the system adding the necessary counterions.'}
    iso = {'type': 'boolean', 'description': '(“False”) Make the box isometric.'}
    positive_ions_number = {'type': 'number', 'description':
                            '(0) Number of additional positive ions to include in the system box.'}

    negative_ions_number = {'type': 'number', 'description':
                            '(0) Number of additional negative ions to include in the system box.'}
    positive_ions_type = {'type': 'string', 'description':
                          '(“Na+”) Type of additional positive ions to include in the system box. Values: Na+,K+.'}
    negative_ions_type = {'type': 'string', 'description':
                          '(“Cl-”) Type of additional negative ions to include in the system box. Values: Cl-.'}

    distance_to_molecule = {'type': 'number', 'description':
                            '(“8.0”) Size for the MD system box -in Angstroms-, defined such as the minimum distance between'
                            'any atom originally present in solute and the edge of the periodic box is given by this distance parameter.'}
    closeness = {'type': 'number', 'description': '(“1.0”) How close, in Å, solvent ATOMs may come to solute ATOMs.'}

    schema = default_schema()
    schema['properties'] = {'forcefield': forcefield, 'water_type': water_type, 'box_type': box_type,
                            'ions_type': ions_type, 'neutralise': neutralise, 'iso': iso,
                            'positive_ions_number': positive_ions_number, 'negative_ions_number': negative_ions_number,
                            'positive_ions_type': positive_ions_type, 'negative_ions_type': negative_ions_type,
                            'distance_to_molecule': distance_to_molecule, 'closeness': closeness}

    return schema


def biobb_leap_add_ions_schema() -> Json:
    """A schema for leap_add_ions options.

    Returns:
        Json: A schema for leap_add_ions options.
    """
    # See https://biobb-amber.readthedocs.io/en/latest/leap.html#module-leap.leap_add_ions
    forcefield_terms = ['protein.ff14SB', 'protein.ff19SB', 'DNA.bsc1', 'DNA.OL15', 'RNA.OL3', 'gaff']
    forcefield = {'type': 'array', 'items': {'type': 'string', 'enum': forcefield_terms},
                  'description': 'Forcefield to be used for the structure generation.'}
    water_type_terms = ['POL3BOX', 'QSPCFWBOX', 'SPCBOX', 'SPCFWBOX', 'TIP3PBOX', 'TIP3PFBOX',
                        'TIP4PBOX', 'TIP4PEWBOX', 'OPCBOX', 'OPC3BOX', 'TIP5PBOX']
    water_type = {'type': 'string', 'enum': water_type_terms,
                  'description': 'Water molecule parameters to be used for the topology'}
    box_type = {'type': 'string', 'enum': ['cubic', 'truncated_octahedron'],
                'description': '(“truncated_octahedron”) Type for the MD system box'}
    ions_type_terms = ['ionsjc_tip3p', 'ionsjc_spce', 'ionsff99_tip3p', 'ions_charmm22', 'ionsjc_tip4pew', 'None']
    ions_type = {'type': 'string', 'enum': ions_type_terms,
                 'description': '(“ionsjc_tip3p”) Ions type.'}
    neutralise = {'type': 'boolean', 'description':
                  '(“False”) Energetically neutralise the system adding the necessary counterions.'}
    ionic_concentration = {'type': 'number', 'description':
                           '(50) Additional ionic concentration to include in the system box. Units in Mol/L.'}
    positive_ions_number = {'type': 'number', 'description':
                            '(0) Number of additional positive ions to include in the system box.'}

    negative_ions_number = {'type': 'number', 'description':
                            '(0) Number of additional negative ions to include in the system box.'}
    positive_ions_type = {'type': 'string', 'description':
                          '(“Na+”) Type of additional positive ions to include in the system box. Values: Na+,K+.'}
    negative_ions_type = {'type': 'string', 'description':
                          '(“Cl-”) Type of additional negative ions to include in the system box. Values: Cl-.'}

    schema = default_schema()
    schema['properties'] = {'forcefield': forcefield, 'water_type': water_type, 'box_type': box_type,
                            'ions_type': ions_type, 'neutralise': neutralise,
                            'ionic_concentration': ionic_concentration,
                            'positive_ions_number': positive_ions_number, 'negative_ions_number': negative_ions_number,
                            'positive_ions_type': positive_ions_type, 'negative_ions_type': negative_ions_type}

    return schema


energy_terms_amber = ['VOLUME', 'TSOLVENT', 'TSOLUTE', 'TEMP', 'PRES',
                      'ETOT', 'ESCF', 'EPTOT', 'EKTOT', 'EKCMT', 'DENSITY']
term_amber = {'type': 'string', 'enum': energy_terms_amber}
terms_amber = {'type': 'array', 'items': term_amber, 'description': '([“ETOT”]) Statistics descriptors.'}
biobb_process_mdout_schema = {**default_schema(), 'properties': {'terms': terms_amber}}


config_schemas = {
    'pdb': biobb_pdb_schema(),
    'pdb2gmx': biobb_pdb2gmx_schema(),
    'editconf': biobb_editconf_schema(),
    'solvate': biobb_solvate_schema(),
    'genion': biobb_genion_schema(),
    'grompp': biobb_grompp_schema(),
    'mdrun': biobb_mdrun_schema(),
    'gmx_energy': biobb_gmx_energy_schema(),
    'gmx_rms': biobb_selection_schema(),
    'gmx_rgyr': biobb_selection_schema(),
    'extract_model': {**default_schema(),
                      'properties': {'models': {'type': 'array', 'items': {'type': 'integer', 'minimum': 0}}}},  # 1?
    'extract_model_pdbqt': {**default_schema(), 'properties': {'model': {'type': 'integer', 'minimum': 0}}},  # 1?
    'str_check_add_hydrogens': biobb_str_check_add_hydrogens_schema(),
    'pdb4amber_run': biobb_pdb4amber_run_schema(),
    'process_mdout': biobb_process_mdout_schema,
    'cpptraj_rms': biobb_cpptraj_rms_schema(),
    'cpptraj_rms_nofit': biobb_cpptraj_rms_nofit_schema(),
    'cpptraj_rgyr': biobb_cpptraj_rgyr_schema(),
    'cpptraj_strip_ambmask': cpptraj_strip_ambmask_schema(),
    'cpptraj_convert': cpptraj_selection_schema(),
    'make_ndx': make_ndx_schema(),
    'genrestr': genrestr_schema(),
    'gmx_trjconv_str': gmx_trjconv_str_schema(),
    'leap_gen_top': leap_gen_top_schema(),
    'process_minout': biobb_process_minout_schema(),
    'leap_solvate': biobb_leap_solvate_schema(),
    'leap_add_ions': biobb_leap_add_ions_schema(),
    'cpptraj_image': cpptraj_selection_schema()
}
