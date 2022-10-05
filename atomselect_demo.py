import mdtraj

def main(selection_string, input_pdb_file):
    traj = mdtraj.load(input_pdb_file)
    print(traj)
    selection_indices = traj.topology.select(selection_string)
    print(selection_indices)
    traj.restrict_atoms(selection_indices)
    traj.save('selection.pdb')


from workflow_types import *

inputs = {'selection_string': string,
          'input_pdb_file': pdbfile}
outputs = {'output_pdb_file': ('selection.pdb', pdbfile)}
