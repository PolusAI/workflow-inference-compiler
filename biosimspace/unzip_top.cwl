#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: unzips a gromacs topology TOP file (and posre.itp).

doc: |-
  unzips a gromacs topology TOP file (and posre.itp).

baseCommand: unzip

inputs:
  input_top_zip_path:
    label: Input zip file
    doc: |-
      Input zip file
      Type: string
      File type: input
      Accepted formats: zip
      Example file: https://github.com/bioexcel/biobb_md/blob/master/biobb_md/test/data/gromacs/genion.zip
    type: File
    format:
    - edam:format_3987
    inputBinding:
      position: 1

outputs:
  output_top_path:
    label: Output topology file
    doc: |-
      Output topology file
      Type: string
      File type: output
      Format: top
      Example file: https://github.com/bioexcel/biobb_md/raw/master/biobb_md/test/data/gromacs/mdrun.top
    type: File
    format: edam:format_3880
    outputBinding:
      glob: p2g.top

# TODO: Add posre.itp

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl