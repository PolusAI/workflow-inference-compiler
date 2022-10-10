# pylint: disable=import-outside-toplevel,no-member
# type: ignore
import distutils.util
import math
from typing import List

def calculate_dG(Kd):
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


def main(input_excel_path, query, min_row, max_row, smiles_column,
         binding_data_column, convert_Kd_dG, output_txt_path):
    import pandas
    import rdkit
    from rdkit import Chem
    from rdkit.Chem import AllChem

    df = pandas.read_excel(input_excel_path, sheet_name=1) # Requires openpyxl

    print(df.shape)
    print(df.columns)
    #print(df)
    print()

    # Determine the categorical column values for query purposes
    #for col in df.columns:
    #    s = list(set(df[col]))
    #    if len(s) < 100: # If more than 100 unique values, column is probably not categorical
    #        print(col, s)

    # For ncats_phenotypic_curated.csv
    # duplicate-classifier ['duplicate', 'unique']
    # Virus ['Dengue Virus', 'Sandfly_Fever', 'HCoV-229E', 'MERS-CoV', 'Yellow Fever Virus', 'Zika Virus', 'RSV', 'Powassan', 'SARS-CoV-2', 'H7N7', 'H1N2', 'HPIV-3', 'West Nile Virus']
    # BAO Label ['cell-based format']
    # Cell_Type [nan, 'Unknown', 'Huh-7', 'Hep-2', 'HEL', 'CEF', 'HG23', 'HMC3', 'HBMEC', 'BHK-21', 'Hepa1-6', 'K-562', 'BHK1', 'J774A.1', 'BHK15', 'MA104', 'BHK', 'LLC-MK2', 'BHK21', 'A549/BHK21', 'MK2', 'JEG3', 'CaCo-2', 'MT-4', 'Hela', 'kidney', 'Vero 76', 'WS1', 'HuH-7', 'Vero', 'HEP-2', 'RAW 264.7', 'Huh-5-2', 'C6/36', 'BHK-D2RepT', 'HuH7', 'MDCK', 'NSC', 'HelaM', 'PBMC', 'A549', 'HELF', 'HEK293', 'BSC-40', 'HAE', 'HFF', 'TREx293', 'MDDC', 'BE(2)-C', 'EAC', 'Vero C1008', 'PEK', 'CEM', 'HEK-293', 'Caco-2', 'BHK-WII', 'HeLa', 'HepG2', 'Vero-76']
    # Standard Type ['EC90', 'Activity', 'CC50', 'EC50', 'TD50', 'Cytotoxicity', 'ID50', 'MIC', 'IC90', 'MCC50', 'Dose', 'Inhibition', 'EC99', 'pIC50', 'MNTD', 'ED50', 'MIC50', 'IC50']
    # Standard Relation [nan, "'~'", "'<'", "'<='", "'>'", "'='", "'>='"]
    # Standard Units [nan, 'uM', '%']
    # Outcome ['Active', 'Inactive', 'Inconclusive']
    # Assay_Type [nan, 'Viral_Replication', 'Unknown', 'Cell_Viability', 'Plaque_Inhibition', 'Focus_Reduction_Assay', 'Proliferation', 'Antigen_Expression', 'Staining_Based', 'Flourescence', 'Viral_Titer', 'Cell_Viability_By_Neutral_Red_Uptake', 'eGFP_Reduction', 'Immunofluorescence', 'Protein_Expression', 'CFI', 'Green_Flourescent_Protein_(eGFP)', 'Viral_Infection', 'Microscopy', 'Immunodetection', 'Replicon_Assay', 'Antigen_Synthesis', 'Viral_RNA_Detection,Plaque_Inhibition,Cell_Viability', 'Cell-based_flavivirus_infection_(CFI)_assay', 'Viral_Yield_Reduction', 'Focus_Forming_Unit_(FFU)_Assay', 'Luciferase', 'Viral_RNA_Detection', 'Luciferase_Reporter_Assay', 'Viral_Entry', 'MTT_Assay', 'RT-PCR', 'Cytopathy', 'Flow_Cytometry', 'Colorimetric', 'Luciferase_Reporter_Gene', 'Cell_Titer', 'Western_Blot', 'Cytotoxicity', 'SDS-PAGE', 'Fluorescence', 'Image-Based', 'Crystal_Violet_Staining_Assay', 'Viral_Reduction_Assay']

    # For ncats_target_based_curated.csv
    # duplicate-type-classifier ['unique', 'duplicate']
    # Virus ['SNV', 'Zika', 'West_Nile', 'RSV', 'SARS-CoV-2', 'H7N7', '229E', 'MERS-CoV', 'HPIV3', 'Dengue']
    # Target Type ['PROTEIN COMPLEX', 'UNCHECKED', 'ORGANISM', 'SINGLE PROTEIN']
    # Target ['Matrix M2-1', 'Phosphoprotein', 'NS2B-NS3 Protease', 'NS5', 'Integrin alpha-V/beta-3', 'not defined', 'Hemagglutinin-neuraminidase', 'Nucleocapsid protein', 'PLpro', 'Spike protein', 'Main Protease (3CLpro, Mpro)', 'Nucleoprotein', 'Fusion glycoprotein F0', 'Neuraminidase', 'Matrix protein 2', 'RDRP']
    # Outcome ['Inactive', 'Active', 'Unclear', 'Inconclusive', 'Undetermined']
    # Standard Type ['IC50', 'EC90', 'Inhibition', 'Kd', 'EC50', 'Activity', 'Ki']
    # Standard Relation [nan, "<'", "<='", ">'", "='", ">='"]
    # Standard Units [nan, 'uM', '%']
    # Comment [nan, 'active', 'Active', 'Not Active', 'Dtt Insensitive', 'Not Determined']

    df = df.query(query)
    print(df)

    # Perform row slicing (if any)
    if int(min_row) != 1 or int(max_row) != -1:
        # We want to convert to zero-based indices and we also want
        # the upper index to be inclusive (i.e. <=) so -1 lower index.
        df = df[(int(min_row) - 1):int(max_row)]
        print(df)

    # Now restrict to the columns we actually care about.
    columns = [smiles_column, binding_data_column]
    df = df[columns].dropna()
    print(df.shape)
    print(df)

    # Generate 2D and/or 3D conformers
    smiles_binding_data: List[str] = []
    convert_Kd_dG = distutils.util.strtobool(convert_Kd_dG)
    for idx, row in enumerate(df.values):

        (smiles, binding_datum) = row
        microMolar = 0.000001 # uM
        binding_datum = binding_datum * microMolar

        if convert_Kd_dG:
            dG = calculate_dG(binding_datum)
            smiles_binding_data.append(f'{smiles} {binding_datum} {dG}')
        else:
            smiles_binding_data.append(f'{smiles} {binding_datum}')

        # See https://www.rdkit.org/docs/GettingStartedInPython.html#working-with-3d-molecules
        mol_2D: rdkit.Chem.rdchem.Mol = Chem.MolFromSmiles(smiles)
        AllChem.Compute2DCoords(mol_2D)

        mol_3D = Chem.AddHs(mol_2D)
        AllChem.EmbedMolecule(mol_3D)
        AllChem.MMFFOptimizeMolecule(mol_3D)

        filename = f'{idx}.sdf' # chemblid is NOT unique!
        writer = Chem.SDWriter(filename)
        #writer = Chem.rdmolfiles.PDBWriter(filename)
        writer.write(mol_3D)
        writer.close()

    with open(output_txt_path, mode='w', encoding='utf-8') as f:
        f.write('\n'.join(smiles_binding_data))

from workflow_types import *

inputs = {'input_excel_path': xlsxfile,
          'query': string,
          'min_row': {**integer, 'default': 1},
          'max_row': {**integer, 'default': -1},
          'smiles_column': string,
          'binding_data_column': string,
          'convert_Kd_dG': string,
          'output_txt_path': string}
outputs = {'output_txt_path': ('$(inputs.output_txt_path)', textfile),
           'output_sdf_path': ('*.sdf', sdffiles)}
