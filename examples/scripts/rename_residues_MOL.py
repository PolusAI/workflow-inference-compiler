import sys

input_mol2_path = sys.argv[1]
output_mol2_path = sys.argv[2]

# See https://proteopedia.org/wiki/index.php/Amino_Acids
# See https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/primary-sequences-and-the-pdb-format
L_aminos_standard = ['ALA', 'ARG', 'ASN', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
                   'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
L_aminos_titrated = ['ASP', 'HID', 'HIE', 'HIP', 'CYX']
L_aminos = L_aminos_standard + L_aminos_titrated + ['PYL', 'SEC', 'UNL']
Nter_L_aminos = ['N' + a for a in L_aminos]
Cter_L_aminos = ['C' + a for a in L_aminos]

# TODO: Consider D-amino acids?

with open(input_mol2_path, mode='r', encoding='utf-8') as f:
    lines = f.readlines()

lines_new = []
for line in lines:
    l = line
    for resname in L_aminos:
        l = l.replace(resname, 'MOL')
    for resname in Nter_L_aminos + Cter_L_aminos:
        l = l.replace(resname, 'MOL ')
    lines_new.append(l)

with open(output_mol2_path, mode='w', encoding='utf-8') as f:
    f.writelines(lines_new)
