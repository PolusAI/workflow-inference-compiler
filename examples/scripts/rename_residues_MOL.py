import sys

input_mol2_path = sys.argv[1]
output_mol2_path = sys.argv[2]

with open(input_mol2_path, mode='r', encoding='utf-8') as f:
    lines = f.readlines()

lines_new = []
index = 7  # mol2 file format residue name column index
for line in lines:
    l = line
    words = l.split()
    if len(words) >= index:  # TODO: and only for ATOM records
        l = l.replace(words[index], 'MOL')

    lines_new.append(l)

with open(output_mol2_path, mode='w', encoding='utf-8') as f:
    f.writelines(lines_new)
