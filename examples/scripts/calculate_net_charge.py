# pylint: disable=no-member
import argparse
from typing import Optional

from rdkit import Chem
from rdkit.Chem import AllChem


def parse_arguments() -> argparse.Namespace:
    """ This function parses the arguments.

    Returns:
        argparse.Namespace: The command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_path', required=True)
    parser.add_argument('--output_log_path', required=True)
    parser.add_argument('--addhydrogens', required=False, default=False)
    args = parser.parse_args()
    return args

def get_net_charge(file_path: str, addhydrogens:bool=False) -> Optional[int]:
    """ Calculate net charge for a given small ligand

    Args:
        file_path (str): The path to the ligand mol2 format structure
        addhydrogens (bool, optional): Flag for adding hydrogens. Defaults to False.

    Returns:
        int: the calculated net charge
    """

    try:
        mol = Chem.MolFromMol2File(file_path, removeHs=False)
    except Exception:
        return None

    if not mol:
        return None

    if addhydrogens:
        mol = Chem.AddHs(mol)
    try:
        AllChem.ComputeGasteigerCharges(mol,
                                        nIter=50,
                                        throwOnParamFailure=True)
    except Exception:
        return None

    num_atoms = mol.GetNumAtoms()
    net_charge = 0.0
    for atom_idx in range(num_atoms):
        atom = mol.GetAtomWithIdx(atom_idx)
        net_charge += float(atom.GetProp("_GasteigerCharge"))
    return round(net_charge)


def main() -> None:
    """ Reads the command line arguments and calculates net charge
    """
    args = parse_arguments()
    net_charge = get_net_charge(args.input_path,
                                    addhydrogens=args.addhydrogens)
    with open(args.output_log_path, mode='w', encoding='utf-8') as wfile:
        wfile.write(f'Calculated net charge: {net_charge}')

if __name__ == '__main__':
    main()
    