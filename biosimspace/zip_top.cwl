#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: zips a gromacs topology TOP file (and posre.itp).

doc: |-
  zips a gromacs topology TOP file (and posre.itp).

baseCommand: zip
arguments: ['system.zip']

inputs:
  input_top_path:
    label: Input topology file
    doc: |-
      Input topology file
      Type: string
      File type: input
      Accepted formats: top
      Example file: https://github.com/bioexcel/biobb_md/raw/master/biobb_md/test/data/gromacs/mdrun.top
    type: File
    format:
    - edam:format_3880
    inputBinding:
      position: 1

# TODO: Add posre.itp

outputs:
  output_top_zip_path:
    label: Output zip file
    doc: |-
      Output zip file
      Type: string
      File type: output
      Format: zip
      Example file: https://github.com/bioexcel/biobb_md/blob/master/biobb_md/test/data/gromacs/genion.zip
    type: File
    format: edam:format_3987
    outputBinding:
      glob: system.zip

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl