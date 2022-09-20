import argparse
import json
import sys
from typing import List, Tuple, Dict, Any

parser = argparse.ArgumentParser(prog='main', description='Parse one or more autodock_vina log files and filter the decoys by docking score.')
parser.add_argument('--input_log_path')
parser.add_argument('--input_log_paths', nargs='+', default=[])
parser.add_argument('--docking_score_cutoff', type=float)
parser.add_argument('--max_num_decoys_per_ligand', type=int)
parser.add_argument('--max_num_decoys_total', type=int)
args = parser.parse_args()

input_log_path = args.input_log_path
input_log_paths = args.input_log_paths
docking_score_cutoff = args.docking_score_cutoff
max_num_decoys_per_ligand = args.max_num_decoys_per_ligand
max_num_decoys_total = args.max_num_decoys_total

if max_num_decoys_per_ligand == -1:
    max_num_decoys_per_ligand = sys.maxsize
if max_num_decoys_total == -1:
    max_num_decoys_total = sys.maxsize

if input_log_path:
    with open(input_log_path, mode='r', encoding='utf-8') as f:
        lines = f.readlines()

if input_log_paths:
    lines_all = []
    for path in input_log_paths:
        with open(path, mode='r', encoding='utf-8') as f:
            lines = f.readlines()
        lines_all.extend(lines)
    lines = lines_all

for line in lines:
    print(line)

# After the initial headers, an autodock vina log file consists of
# 1 or more blocks of the following form:
"""
Performing docking (random seed: -193947260) ... 
0%   10   20   30   40   50   60   70   80   90   100%
|----|----|----|----|----|----|----|----|----|----|
***************************************************

mode |   affinity | dist from best mode
     | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
   1       -5.773          0          0
   2       -5.577      4.821      7.207
   3        -5.46      4.053      6.446
   4        -5.44      1.344      3.349
   5       -5.431      3.576       5.58
   6       -5.412      3.012      5.868
   7       -5.405      1.914      4.145
   8       -5.403      3.209      6.271
   9       -5.392      2.556      4.714
"""

scores: List[float] = []
scores_all: List[List[float]] = []
parsing = False
for line in lines:
    if line.startswith('-----+------------+----------+----------'):
        scores = []
        parsing = True
        continue

    if parsing:
        try:
            strs = line.split()
            mode_idx = int(strs[0])
            floats = [float(x) for x in strs[1:]]
            score = floats[0]
            if score < docking_score_cutoff and len(scores) < max_num_decoys_per_ligand:
                scores.append(score)
        except Exception as e:
            scores_all.append(scores)
            parsing = False

if parsing: # When we reach end of file, we need to append one last time.
    scores_all.append(scores)

# Now we need to globally sort (by the docking score) and apply the max total, while preserving the 2D structure.
indexed_scores: List[Tuple[float, Tuple[int, int]]] = []
for mol_idx, scores in enumerate(scores_all):
    for mode_idx, score in enumerate(scores):
        indexed_scores.append((score, (mol_idx, mode_idx)))

indexed_scores.sort(key=lambda x: x[0]) # Sort by the docking scores
indexed_scores = indexed_scores[:max_num_decoys_total] # Truncate upto max total
#indexed_scores.sort(key=lambda x: x[1]) # Sort by the index tuple, which uses dictionary order.

indices_all: List[str] = []
for (score, (mol_idx, mode_idx)) in indexed_scores:
    indices_all.append(str(score) + " " + str(mol_idx) + " " + str(mode_idx))

with open('indices.txt', mode='w', encoding='utf-8') as f:
    f.write('\n'.join(indices_all))

# TODO: Since outputEval doesn't seem to work for outputting 2D arrays,
# Try to specify the outputs by writing cwl.output.json
# See https://cwl.discourse.group/t/cwl-output-json/579
#outputs: Dict[str, Any] = {}
#outputs['output_batch_pdbqt_path'] = [[]]

#with open('cwl.output.json', mode='w', encoding='utf-8') as f:
#    f.write(json.dumps(outputs))
