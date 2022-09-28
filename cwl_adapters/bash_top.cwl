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

  input_top_path:
    type: File
    format:
    - edam:format_3880
    inputBinding:
      position: 2

  output_top_path:
    type: string
    format:
    - edam:format_2330 # 'Textual format'
#    inputBinding:
#      position: 3
    default: system.top

outputs:
  output_top_path:
    type: File
    format: edam:format_3880
    streamable: true
    outputBinding:
      glob: $(inputs.output_top_path)

stdout: $(inputs.output_top_path)

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
