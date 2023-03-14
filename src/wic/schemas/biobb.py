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
    # type "MainChain", IntelliSense correctly suggests "MainChain", "MainChain+Cb",
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

    schema = default_schema()
    schema['properties'] = {'distance_to_molecule': distance_to_molecule,
                            'box_type': box_type, 'center_molecule': center_molecule,
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
    schema['properties'] = {'water_type': water_type, 'forcefield': forcefield,
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


energy_terms_amber = ['VOLUME', 'TSOLVENT', 'TSOLUTE', 'TEMP', 'PRES',
                      'ETOT', 'ESCF', 'EPTOT', 'EKTOT', 'EKCMT', 'DENSITY']
term_amber = {'type': 'string', 'enum': energy_terms_amber}
terms_amber = {'type': 'array', 'items': term_amber, 'description': '([“ETOT”]) Statistics descriptors.'}
biobb_process_mdout_schema = {**default_schema(), 'properties': {'terms': terms_amber}}


config_schemas = {
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
}
