#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: Concatenate multiple sdf files into a single sdf file

doc: |
  Concatenate multiple sdf files into a single sdf file

baseCommand: cat

inputs: 
  input_sdfs_path:
    label: Input sdf files
    doc: |-
      Input sdf files
      Type: File[]
    type: File[]
    format:
    - edam:format_3814
    inputBinding:
      position: 1

  output_sdf_path:
    type: string
    format:
    - edam:format_3814
    default: system.sdf

outputs:
  output_sdf_path:
    type: File
    format: edam:format_3814
    streamable: true
    outputBinding:
      glob: $(inputs.output_sdf_path)

stdout: $(inputs.output_sdf_path)

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl