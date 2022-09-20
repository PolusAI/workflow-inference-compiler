#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: This class is a wrapper of the Open Babel tool.

doc: |-
  Small molecule format conversion for structures or trajectories. Open Babel is a chemical toolbox designed to speak the many languages of chemical data. It's an open, collaborative project allowing anyone to search, convert, analyze, or store data from molecular modeling, chemistry, solid-state materials, biochemistry, or related areas. Visit the official page.

baseCommand: obabel
# "-r	Remove all but the largest contiguous fragment (strip salts)"

# NOTE: The version of openbabel in this container is old; This may or may not cause problems.
# See https://github.com/openbabel/openbabel/issues/2435
hints:
  DockerRequirement:
    #dockerPull: jakefennick/scripts
    dockerPull: quay.io/biocontainers/biobb_chemistry:3.7.0--pyhdfd78af_0

inputs:
  first_molecule:
    label: Index of the first molecule (1-based)
    doc: |-
      Input Index of the first molecule (1-based)
      Type: string
    type: string?
    format:
    - edam:format_2330 # 'Textual format'
    inputBinding:
      prefix: -f
    default: 1

  last_molecule:
    label: Index of the last molecule (1-based)
    doc: |-
      Input Index of the last molecule (1-based)
      Type: string
    type: string?
    format:
    - edam:format_2330 # 'Textual format'
    inputBinding:
      prefix: -l
    default: 1

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
    label: Path to the output file
    doc: |-
      Path to the output file
      Type: string
      File type: output
      Accepted formats: sdf
    type: string
    format:
    - edam:format_3814 # sdf
    inputBinding:
      position: 2
      prefix: -O
    default: system.sdf

  arg1:
    label: Additional arguments
    doc: |-
      Additional arguments
      Type: string
    type: string?
    format:
    - edam:format_2330 # 'Textual format'
    inputBinding:
      position: 3
    default: ""

  arg2:
    label: Additional arguments
    doc: |-
      Additional arguments
      Type: string
    type: string?
    format:
    - edam:format_2330 # 'Textual format'
    inputBinding:
      position: 4
    default: ""

  arg3:
    label: Additional arguments
    doc: |-
      Additional arguments
      Type: string
    type: string?
    format:
    - edam:format_2330 # 'Textual format'
    inputBinding:
      position: 5
    default: ""

  arg4:
    label: Additional arguments
    doc: |-
      Additional arguments
      Type: string
    type: string?
    format:
    - edam:format_2330 # 'Textual format'
    inputBinding:
      position: 6
    default: ""

  arg5:
    label: Additional arguments
    doc: |-
      Additional arguments
      Type: string
    type: string?
    format:
    - edam:format_2330 # 'Textual format'
    inputBinding:
      position: 7
    default: ""

outputs:
  output_sdf_path:
    label: Path to the output file
    doc: |-
      Path to the output file
    type: File
    outputBinding:
      glob: $(inputs.output_sdf_path)
    format: edam:format_3814 # sdf

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
