#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: Convert topology and coordinate files from amber to gromacs format using BioSimSpace.

doc: |-
  Convert topology and coordinate files from amber to gromacs format using BioSimSpace.

baseCommand: python3
arguments: [$(inputs.script), '--output_top_format', 'GroTop', '--output_crd_format', 'Gro87']

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
    label: Input topology file (AMBER ParmTop)
    doc: |-
      Input topology file (AMBER ParmTop)
      Type: string
      File type: input
      Accepted formats: top, parmtop, prmtop
      Example file: https://github.com/bioexcel/biobb_amber/raw/master/biobb_amber/test/data/sander/cln025.prmtop
    type: File
    format:
    - edam:format_3881
    - edam:format_3881
    - edam:format_3881
    inputBinding:
      position: 1
      prefix: --input_top_path

  input_crd_path:
    label: Input coordinates file (AMBER crd)
    doc: |-
      Input coordinates file (AMBER crd)
      Type: string
      File type: input
      Accepted formats: crd, mdcrd, inpcrd, netcdf, nc, ncrst
      Example file: https://github.com/bioexcel/biobb_amber/raw/master/biobb_amber/test/data/sander/cln025.inpcrd
    type: File
    format:
    - edam:format_3878
    - edam:format_3878
    - edam:format_3878
    - edam:format_3650
    - edam:format_3650
    - edam:format_3886
    inputBinding:
      position: 2
      prefix: --input_crd_path

outputs:
  output_top_path:
    label: Output topology file (GROMACS top)
    doc: |-
      Output topology file (GROMACS top)
      Type: string
      File type: output
      Format: top
      Example file: https://github.com/bioexcel/biobb_md/raw/master/biobb_md/test/data/gromacs/mdrun.top
    type: File
    format: edam:format_3880
    outputBinding:
      glob: top.top  # BioSimSpace automatically renames *.GroTop to *.top

  output_crd_path:
    label: Output coordinates file (GROMACS gro)
    doc: |-
      Output coordinates file (GROMACS gro)
      Type: string
      File type: output
      Format: gro
      Example file: https://github.com/bioexcel/biobb_md/raw/master/biobb_md/test/reference/gromacs/ref_mdrun.gro
    type: File
    format: edam:format_2033
    outputBinding:
      glob: top.gro  # BioSimSpace automatically renames *.Gro87 to *.gro

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl