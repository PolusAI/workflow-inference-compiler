#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: This class is a wrapper of Acpype tool for generation of topologies for GROMACS.

doc: |-
  Generation of topologies for GROMACS. Acpype is a tool based in Python to use Antechamber to generate topologies for chemical compounds and to interface with others python applications like CCPN or ARIA. Visit the official page.

baseCommand: acpype

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
      Accepted formats: pdb, mdl, mol2
      Example file: https://github.com/bioexcel/biobb_chemistry/raw/master/biobb_chemistry/test/data/acpype/acpype.params.mol2
    type: File
    format:
    - edam:format_1476
    - edam:format_3815
    - edam:format_3816
    inputBinding:
      position: 1
      prefix: --input

  output_path_gro:
    label: Path to the GRO output file
    doc: |-
      Path to the GRO output file
      Type: string
      File type: output
      Accepted formats: gro
      Example file: https://github.com/bioexcel/biobb_chemistry/raw/master/biobb_chemistry/test/reference/acpype/ref_acpype.gmx.gro
    type: string
    format:
    - edam:format_2033
    #inputBinding:
    #  position: 2
    #  prefix: --output_path_gro
    default: system.gro

  output_path_itp:
    label: Path to the ITP output file
    doc: |-
      Path to the ITP output file
      Type: string
      File type: output
      Accepted formats: itp
      Example file: https://github.com/bioexcel/biobb_chemistry/raw/master/biobb_chemistry/test/reference/acpype/ref_acpype.gmx.itp
    type: string
    format:
    - edam:format_3883
    #inputBinding:
    #  position: 3
    #  prefix: --output_path_itp
    default: system.itp

  output_path_top:
    label: Path to the TOP output file
    doc: |-
      Path to the TOP output file
      Type: string
      File type: output
      Accepted formats: top
      Example file: https://github.com/bioexcel/biobb_chemistry/raw/master/biobb_chemistry/test/reference/acpype/ref_acpype.gmx.top
    type: string
    format:
    - edam:format_3880
    #inputBinding:
    #  position: 4
    #  prefix: --output_path_top
    default: system.top

  base_name:
    label: Prefix for the output filenames
    doc: |-
      Prefix for the output filenames
      Type: string
    type: string
    format:
    - edam:format_2330
    inputBinding:
      position: 2
      prefix: --basename
    default: ligand

  charge_method:
    label: Method to determine the atomic partial charges
    doc: |-
      Method to determine the atomic partial charges
      Type: string
    type: string
    format:
    - edam:format_2330
    inputBinding:
      position: 3
      prefix: --charge_method
    default: bcc # gas, bcc (default), user (user's charges in mol2 file)

outputs:
  output_path_gro:
    label: Path to the GRO output file
    doc: |-
      Path to the GRO output file
    type: File
    outputBinding:
      glob: "*/$(inputs.base_name)_GMX.gro"
    format: edam:format_2033

  output_path_itp:
    label: Path to the ITP output file
    doc: |-
      Path to the ITP output file
    type: File
    outputBinding:
      glob: "*/$(inputs.base_name)_GMX.itp"
    format: edam:format_3883

  output_path_top:
    label: Path to the TOP output file
    doc: |-
      Path to the TOP output file
    type: File
    outputBinding:
      glob: "*/$(inputs.base_name)_GMX.top"
    format: edam:format_3880

# See https://rabix.io/cwl-patterns.html
#  output_all:
#    type:
#      type: array
#      items: [File]
#    outputBinding: {glob: "*/*"}
#    format: edam:format_2330

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
