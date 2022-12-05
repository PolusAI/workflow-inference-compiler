# pylint: disable=import-outside-toplevel,no-member
# type: ignore
from collections import defaultdict
import distutils.util
import math
import os.path as osp
import re
import subprocess
import argparse

import pandas as pd


def parse_arguments() -> argparse.Namespace:
    """ This function parses the arguments.

    Returns:
        argparse.Namespace: The command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--index_file_name', required=True)
    parser.add_argument('--base_dir', required=True)
    parser.add_argument('--query', required=False, default=False)
    parser.add_argument('--output_txt_path', required=True)
    parser.add_argument('--min_row', required=False, default=1)
    parser.add_argument('--max_row', required=False, default=-1)
    parser.add_argument('--convert_Kd_dG', required=False, default="False")

    args = parser.parse_args()
    return args

def calculate_dG(Kd: float) -> float:
    """ Calculates binding free energy from Kd

    Args:
        Kd (float): The binding affinity of the protein-ligand complex

    Returns:
        float: The binding free energy
    """
    # Calculate the binding free energy from Kd so we can make the correlation plots.
    # See https://en.wikipedia.org/wiki/Binding_constant
    ideal_gas_constant = 8.31446261815324 # J/(Mol*K)
    kcal_per_joule = 4184
    # NOTE: Unfortunately, the temperature at which experimental Kd binding data was taken
    # is often not recorded. Thus, we are forced to guess. The two standard guesses are
    # physiological body temperature (310K) or room temperature (298K).
    temperature = 298
    RT = (ideal_gas_constant / kcal_per_joule) * temperature
    # NOTE: For performance, simulations are often done in a very small unit cell, and
    # thus at a very high concentration. The size of the unit cell bounds the volume.
    # For shorter simulations where the ligand has not explored the entire box, it may
    # be less. See the Yank paper for a method of calculating the correct volumes.
    standard_concentration = 1 # Units of mol / L, but see comment above.
    dG = RT * math.log(Kd / standard_concentration)
    return dG


def read_index_file(index_file_path: str) -> pd.DataFrame:
    """ Reads the PDBbind index file and extracts binding data

    Args:
        index_file_path (str): The path to the index file

    Returns:
        pd.DataFrame: The Kd data
    """
    data = defaultdict(list)
    # The file format
    # PDB code, resolution, release year, -logKd/Ki, Kd/Ki, reference, ligand name
    unit_conv = {'uM': 1,
                 'mM': 1000.0,
                 'nM': 0.001,
                 'pM': 0.000001}

    with open(index_file_path, mode='r', encoding='utf-8') as rfile:
        lines = [line for line in rfile.readlines() if line[0] != '#' and 'Kd=' in line]
        for line in lines:
            words = line.split()
            data['PDB_code'].append(words[0])
            data['resolution'].append(words[1])
            data['release_year'].append(words[2])

            # Kd conversion to micro molar
            unit = re.split(r"=[-+]?(?:\d*\.\d+|\d+)", words[4])[1]
            standard_type = re.split(r"=[-+]?(?:\d*\.\d+|\d+)", words[4])[0]
            kd = float(re.findall(r"[-+]?(?:\d*\.\d+|\d+)",words[4])[0])
            data['Kd_Ki'].append(standard_type)
            data['value'].append(kd * unit_conv[unit])
            data['ligand_name'].append(re.findall(r'\((.*?)\)', words[7])[0])

    return pd.DataFrame.from_dict(data)


def load_data(index_file_name: str, base_dir: str, query:str, output_txt_path:str,
              min_row: int=1, max_row: int=-1, convert_Kd_dG: bool=False) -> None:
    """ Filters Kd data beased on a query

    Args:
        index_file_name (str): The PDBbind index file name
        base_dir (str): The base directry of the dataset
        query (str): The Query to perform
        output_txt_path (str): The ouput text file
        min_row (int, optional): min index of rows. Defaults to 1.
        max_row (int, optional): max index of rows. Defaults to -1.
        convert_Kd_dG (bool, optional): If this set to True, The dG will be calculated. Defaults to False.
    """

    index_file_path = osp.join(base_dir, 'index', index_file_name)
    df = read_index_file(index_file_path)
    print(df.shape)
    print(df.columns)
    # perform query
    df = df.query(query)

    # Perform row slicing (if any)
    if int(min_row) != 1 or int(max_row) != -1:
        # We want to convert to zero-based indices and we also want
        # the upper index to be inclusive (i.e. <=) so -1 lower index.
        df = df[(int(min_row) - 1):int(max_row)]

    # Calculate dG
    convert_Kd_dG = distutils.util.strtobool(convert_Kd_dG)
    binding_data = df[['PDB_code', 'value', 'Kd_Ki']]
    if convert_Kd_dG:
        microMolar = 0.000001  # uM
        dG_data = [calculate_dG(value * microMolar) for value in binding_data['value']]
        binding_data.insert(2, 'dG', dG_data)

    with open(output_txt_path, mode='w', encoding='utf-8') as f:
        dfAsString = binding_data.to_string(header=False, index=False)
        f.write(dfAsString)


    # copy pdb and sdf files
    for _, row in binding_data.iterrows():
        source_pdb_path = osp.join(base_dir,
                                   row['PDB_code'],
                                   f'{row["PDB_code"]}_protein.pdb')
        dist_pdb_path = f'{row["PDB_code"]}_protein.pdb'
        subprocess.run(["cp", f"{source_pdb_path}", f"{dist_pdb_path}"])
        source_sdf_path = osp.join(base_dir,
                                    row['PDB_code'],
                                   f'{row["PDB_code"]}_ligand.sdf')

        dist_sdf_path = f'{row["PDB_code"]}_ligand.sdf'
        subprocess.run(["cp", f"{source_sdf_path}", f"{dist_sdf_path}"])

def main() -> None:
    """ Reads the command line arguments
    """
    args = parse_arguments()
    load_data(args.index_file_name, args.base_dir, args.query, args.output_txt_path,
              min_row=args.min_row, max_row=args.max_row, convert_Kd_dG=args.convert_Kd_dG)

if __name__ == '__main__':
    main()