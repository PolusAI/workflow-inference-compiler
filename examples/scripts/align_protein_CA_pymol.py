import subprocess as sub
import sys

import pymol

# The `polymer` selection keyword incorrectly assigns the C-terminal amino
# acid residue to the ligand. Alternatively, we can select atoms based on their
# indices. Thus, these two files are only needed to get the number of atoms.
receptorpdb = sys.argv[1]
ligandpdb = sys.argv[2]

geniongro = sys.argv[3]
prodtrr = sys.argv[4]
savepdb = sys.argv[5]

pymol.cmd.load(receptorpdb, 'receptor')
pymol.cmd.load(ligandpdb, 'ligand')

R = pymol.cmd.select('R', 'receptor')
# This assumes the receptor atoms are first in the file.
indices_receptor = ' or '.join([f'index {i}' for i in range(1, 1 + int(R))])

L = pymol.cmd.select('L', 'ligand')
# This assumes the ligand atoms are appended to the file after the receptor.
indices_ligand = ' or '.join([f'index {i}' for i in range(1 + int(R), 1 + int(R) + int(L))])


pymol.cmd.load(geniongro, 'genion')
G = pymol.cmd.select('G', f'genion and name CA and ({indices_receptor})')

# NOTE: We do NOT want to load prod.gro here because the coordinates only get
# written after the last timestep, so we can't use it for realtime analysis.
# However, pymol just needs to extract the topology information, so we can use
# genion.gro and overwrite the coordinates from the trajectory file prod.trr
pymol.cmd.load(geniongro, 'prod')
interval = 1 # Only load every nth frame.

# NOTE: Pymol performs an alignment using an average across all frames in the trajectory;
# See https://sourceforge.net/p/pymol/mailman/message/29514454/
# It does NOT perform a separate alignment for each frame!
# (That's what we really want because e.g. gromacs cannot remove the angular drift
# when using a periodic simulation box. It can still remove the linear drift, however.
# This is probably fine for very short simulations, but eventually the angular drift may
# drift may inflate the rmsd of the ligand w.r.t. the docked conformation.
# This is bad because we really want to avoid false negatives.)

ALIGN_SEPARATELY = False
if ALIGN_SEPARATELY:
    # We can process each frame individually and concatenate the pdb files together.
    # Unfortunately, re-seeking into the trajectory file causes this loop to be O(n^2)
    sub.run(['rm', savepdb], check=False) # Delete savepdb (if it exists) due to >> below.

    # First load all of the states so we can call count_states.
    pymol.cmd.load(geniongro, 'num_states')
    pymol.cmd.load_traj(prodtrr, 'num_states')
    num_states = pymol.cmd.count_states('num_states')
    #print('num_states', num_states)

    #states = range(1, 1 + num_states)
    #interval = max(1, round(num_states/1000)) # More than 1000 = SLOW
    states = [interval*i for i in range(1, 1 + int(num_states/interval))]

    pymol.cmd.load(geniongro, 'aligned_sep')
    for idx, state in enumerate(states):
        # Now re-load each state individually.
        pymol.cmd.load_traj(prodtrr, 'prod', state=1, start=state, stop=state, max=1)

        P = pymol.cmd.select('P', f'prod and name CA and ({indices_receptor})')
        F = pymol.cmd.pair_fit('P', 'G')
        #print(F)

        pymol.cmd.select('N', f'prod and ({indices_receptor} or {indices_ligand})')
        pymol.cmd.save(f'{state}.pdb', 'N', -1) # 0 == all states (default -1 == current state only)

        # For interactive visualization / comparison, load the aligned results.
        #pymol.cmd.select('T', 'prod')
        pymol.cmd.save('temp.pdb', 'prod', -1) # 0 == all states (default -1 == current state only)
        pymol.cmd.load_traj('temp.pdb', 'aligned_sep', state=1+idx) # 0 = append
        sub.run(['rm', 'temp.pdb'], check=False)

        # Now concatenate the individual pdb files into savepdb
        # TODO: For security, avoid using shell=True
        sub.run(f'cat {state}.pdb >> {savepdb}', shell=True, check=True)
        sub.run(['rm', f'{state}.pdb'], check=False)

    # For interactive visualization / comparison, load the aligned results.
    #pymol.cmd.load(savepdb, 'aligned_sep')
    #pymol.cmd.load_traj(savepdb, 'aligned_sep')

ALIGN_AVERAGED = True
if ALIGN_AVERAGED:
    # Load all of the states.
    pymol.cmd.load_traj(prodtrr, 'prod', interval=interval)

    P = pymol.cmd.select('P', f'prod and name CA and ({indices_receptor})')
    F = pymol.cmd.pair_fit('P', 'G')
    #print(F)

    pymol.cmd.select('N', f'prod and ({indices_receptor} or {indices_ligand})')
    pymol.cmd.save(savepdb, 'N', 0) # 0 == all states (default -1 == current state only)

MDTRAJ = True
if MDTRAJ:
    import MDAnalysis
    from MDAnalysis.analysis import align
    from MDAnalysis.coordinates import TRR

    gro = MDAnalysis.Universe(topology=geniongro, coordinates=geniongro)
    print(gro)
    # NOTE: Setting coordinates= in the constructor only loads the first frame!
    trr = MDAnalysis.Universe(topology=geniongro)
    trr.trajectory = TRR.TRRReader(prodtrr) # This loads all frames.
    print(trr)

    # Don't forget to call .run()! Otherwise, it will silently do nothing
    # (except write an empty file).
    aligntraj = align.AlignTraj(trr, gro, select='protein and name CA', filename='aligned.trr').run()
    print(aligntraj.frames)

    pymol.cmd.load(geniongro, 'aligned_mda')
    pymol.cmd.load_traj('aligned.trr', 'aligned_mda')
