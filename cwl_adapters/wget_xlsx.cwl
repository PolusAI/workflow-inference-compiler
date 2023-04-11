#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: Download an xlsx file

doc: |
  Download an xlsx file

baseCommand: wget

requirements:
  InlineJavascriptRequirement: {}

inputs:
  url:
    label: url
    doc: |-
      url
      Type: string
    type: string
    inputBinding:
      position: 1

  output_xlsx_path:
    type: string
    format:
    - edam:format_2330
    inputBinding:
      position: 2
      prefix: -O
    default: system.xlsx

outputs:
  output_xlsx_path:
    type: File
    format: edam:format_3620
    streamable: true
    outputBinding:
      glob: $(inputs.output_xlsx_path)

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl