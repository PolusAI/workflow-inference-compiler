#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: Run a python3 script

doc: |-
  Run a python3 script

baseCommand: python3

hints:
  DockerRequirement:
    dockerPull: jakefennick/scripts

inputs:
  script:
    type: string
    inputBinding:
      position: 1

  input_mol2_path:
    type: File
    format:
    - edam:format_3816 # mol2
    inputBinding:
      position: 2

  output_mol2_path:
    label: Path to the output file
    doc: |-
      Path to the output file
    type: string
    format:
    - edam:format_3816 # mol2
    inputBinding:
      position: 3
    default: system.mol2

outputs:
  output_mol2_path:
    label: Path to the output file
    doc: |-
      Path to the output file
    type: File
    format: edam:format_3816 # mol2
    outputBinding:
      glob: $(inputs.output_mol2_path)

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl