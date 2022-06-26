import sys

import MDAnalysis

input_mol2 = sys.argv[1]
output_pdb = sys.argv[2]

u = MDAnalysis.Universe(input_mol2)

writer = MDAnalysis.coordinates.PDBQT.PDBQTWriter(output_pdb) # ligand.pdbqt'
writer.write(u)