#!/usr/bin/env cwl-runner
cwlVersion: v1.2

class: CommandLineTool

label: Uses openbabel to add hydrogens and minimize a small molecule, search for the lowest energy conformer, then minimize again.

doc: |
  Uses openbabel to add hydrogens and minimize a small molecule, search for the lowest energy conformer, then minimize again.

baseCommand: obgen

# NOTE: The version of openbabel in this container is old; This may or may not cause problems.
# See https://github.com/openbabel/openbabel/issues/2435
hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/biobb_chemistry:3.7.0--pyhdfd78af_0

inputs: 
  input_path:
    label: Path to the input file
    doc: |-
      Path to the input file
      Type: string
      File type: input
      Accepted formats: dat, ent, fa, fasta, gro, inp, log, mcif, mdl, mmcif, mol, mol2, pdb, pdbqt, png, sdf, smi, smiles, txt, xml, xtc
      Example file: https://github.com/bioexcel/biobb_chemistry/raw/master/biobb_chemistry/test/data/babel/babel.smi
    type: File
    format:
    - edam:format_1637
    - edam:format_1476
    - edam:format_1929
    - edam:format_1929
    - edam:format_2033
    - edam:format_3878
    - edam:format_2030
    - edam:format_1477
    - edam:format_3815
    - edam:format_1477
    - edam:format_3815
    - edam:format_3816
    - edam:format_1476
    - edam:format_1476
    - edam:format_3603
    - edam:format_3814
    - edam:format_1196
    - edam:format_1196
    - edam:format_2033
    - edam:format_2332
    - edam:format_3875
    inputBinding:
      position: 1

  output_sdf_path:
    type: string
    format:
    - edam:format_3814 # sdf

outputs:
  output_sdf_path:
    type: File
    format: edam:format_3814 # sdf
    streamable: true
    outputBinding:
      glob: $(inputs.output_sdf_path)

stdout: $(inputs.output_sdf_path)

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
