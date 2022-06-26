#!/usr/bin/env cwl-runner
cwlVersion: v1.2

class: CommandLineTool

label: Run a Bash script

doc: |
  Run a Bash script

baseCommand: bash

inputs: 
  script:
    type: File
    format:
    - edam:format_2330 # 'Textual format'
    inputBinding:
      position: 1

  input_pdb_path:
    type: File
    format:
    - edam:format_1476
    inputBinding:
      position: 2

  output_pdb_path:
    type: string
    format:
    - edam:format_1476
#    inputBinding:
#      position: 3

outputs:
  output_pdb_path:
    type: File
    format: edam:format_1476
    streamable: true
    outputBinding:
      glob: $(inputs.output_pdb_path)

stdout: $(inputs.output_pdb_path)

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
