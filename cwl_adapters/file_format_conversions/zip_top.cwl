#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: zips a gromacs topology TOP file (and/or itp include file).

doc: |-
  zips a gromacs topology TOP file (and/or itp include file).

baseCommand: zip
arguments: ["-j", $(inputs.output_top_zip_path)] # junk (don't record) directory names

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

  input_itp_path:
    label: Input topology include file
    doc: |-
      Input topology include file
      Type: string
      File type: input
      Accepted formats: top
      Example file: https://github.com/bioexcel/biobb_md/raw/master/biobb_md/test/data/gromacs/mdrun.top
    type: File
    format:
    - edam:format_3883
    inputBinding:
      position: 2

  output_top_zip_path:
    label: Output zip file
    type: string
    format: edam:format_2330 # 'textual format'
    default: system.zip

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
      glob: $(inputs.output_top_zip_path)

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl