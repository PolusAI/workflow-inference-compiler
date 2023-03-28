from workflow_types import *
# NOTE: No other top-level imports supported


def main(selection_string, input_pdb_path, output_pdb_path):   # type: ignore[no-untyped-def]
    import mdtraj
    traj = mdtraj.load(input_pdb_path)
    print(traj)
    selection_indices = traj.topology.select(selection_string)
    print(selection_indices)
    traj.restrict_atoms(selection_indices)
    traj.save(output_pdb_path)


inputs = {'selection_string': string,
          'input_pdb_path': pdbfile,
          'output_pdb_path': {**string, 'default': 'selection.pdb'}}
outputs = {'output_pdb_path': ('$(inputs.output_pdb_path)', pdbfile)}
