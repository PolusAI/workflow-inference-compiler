#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: Convert topology and coordinate files from gromacs to amber format using BioSimSpace.

doc: |-
  Convert topology and coordinate files from gromacs to amber format using BioSimSpace.

baseCommand: python3
arguments: [$(inputs.script), '--output_top_format', 'PRM7', '--output_crd_format', 'RST7']

hints:
  DockerRequirement:
    dockerPull: jakefennick/biosimspace

inputs:
  script:
    type: string
    inputBinding:
      position: 1
    default: /conversion.py # NOTE: Initial / required

  input_top_path:
    label: Path to the portable binary run input file TOP
    doc: |-
      Path to the portable binary run input file TOP
      Type: string
      File type: input
      Accepted formats: top
      Example file: https://github.com/bioexcel/biobb_md/raw/master/biobb_md/test/data/gromacs/mdrun.top
    type: File
    format:
    - edam:format_3880
    inputBinding:
      position: 1
      prefix: --input_top_path

  input_crd_path:
    label: Path to the input GROMACS structure GRO file
    doc: |-
      Path to the input GROMACS structure GRO file
      Type: string
      File type: input
      Accepted formats: gro
      Example file: https://github.com/bioexcel/biobb_md/raw/master/biobb_md/test/reference/gromacs/ref_mdrun.gro
    type: File
    format:
    - edam:format_2033
    inputBinding:
      position: 2
      prefix: --input_crd_path

outputs:
  output_top_path:
    label: Output topology file (AMBER ParmTop)
    doc: |-
      Output topology file (AMBER ParmTop)
      Type: string
      File type: output
      Format: prmtop
      Example file: https://github.com/bioexcel/biobb_amber/raw/master/biobb_amber/test/data/sander/cln025.prmtop
    type: File
    format: edam:format_3881
    outputBinding:
      glob: top.prmtop

  output_crd_path:
    label: Output coordinates file (AMBER crd)
    doc: |-
      Output coordinates file (AMBER crd)
      Type: string
      File type: output
      Format: crd
      Example file: https://github.com/bioexcel/biobb_amber/raw/master/biobb_amber/test/data/sander/cln025.inpcrd
    type: File
    format: edam:format_3878
    outputBinding:
      glob: top.inpcrd

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl