#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: Run a Bash script on two pdb files

doc: |
  Run a Bash script on two pdb files

baseCommand: bash

hints:
  DockerRequirement:
    dockerPull: jakefennick/scripts

inputs: 
  script:
    type: string
    inputBinding:
      position: 1

  input_pdb1_path:
    type: File
    format:
    - edam:format_1476
    inputBinding:
      position: 2

  input_pdb2_path:
    type: File
    format:
    - edam:format_1476
    inputBinding:
      position: 3

  output_path:
    type: string
    default: system.txt

outputs:
  output_path:
    type: File
    streamable: true
    outputBinding:
      glob: $(inputs.output_path)

stdout: $(inputs.output_path)

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
