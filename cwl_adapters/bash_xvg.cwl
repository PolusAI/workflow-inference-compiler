#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: Run a Bash script

doc: |
  Run a Bash script

baseCommand: bash

hints:
  DockerRequirement:
    dockerPull: jakefennick/scripts

inputs: 
  script:
    type: string
    inputBinding:
      position: 1

  input_xvg1_path:
    type: File
    format:
    - edam:format_2030
    inputBinding:
      position: 2

  input_xvg2_path:
    type: File
    format:
    - edam:format_2030
    inputBinding:
      position: 2

  output_xvg_path:
    type: string
    format:
    - edam:format_2030
#    inputBinding:
#      position: 3

outputs:
  output_xvg_path:
    type: File
    format: edam:format_2030
    streamable: true
    outputBinding:
      glob: $(inputs.output_xvg_path)

stdout: $(inputs.output_xvg_path)

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
