
#!/usr/bin/env cwl-runner
cwlVersion: v1.2

class: CommandLineTool

label: Calls .toString() to convert any object into a string.

doc: |-
  Calls .toString() to convert any object into a string.

baseCommand: "true"

requirements:
  InlineJavascriptRequirement: {}

inputs:
  input:
    label: The input data
    doc: |-
      The input data
    type: Any
    format:
    - edam:format_2330

outputs:
  output:
    label: The output string
    doc: |-
      The output string
    type: string
    #format: edam:format_2330 # "'str' object does not support item assignment"
    outputBinding:
      outputEval: $(inputs.input.toString())

#stdout: stdout

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
