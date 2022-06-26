import sys

import MDAnalysis
from MDAnalysis.analysis import align
from MDAnalysis.coordinates import TRR

input_gro_path = sys.argv[1]
input_trr_path = sys.argv[2]
output_trr_path = sys.argv[3]

gro = MDAnalysis.Universe(topology=input_gro_path, coordinates=input_gro_path)
print(gro)
# NOTE: Setting coordinates= in the constructor only loads the first frame!
trr = MDAnalysis.Universe(topology=input_gro_path)
trr.trajectory = TRR.TRRReader(input_trr_path) # This loads all frames.
print(trr)

# Don't forget to call .run()! Otherwise, it will silently do nothing
# (except write an empty file).
aligntraj = align.AlignTraj(trr, gro, select='protein and name CA', filename=output_trr_path).run()
print(aligntraj.frames)
