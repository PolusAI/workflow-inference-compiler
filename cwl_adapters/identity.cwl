
#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: The identity operation. Useful for specifying data as a literal array, for scattering in a later step.

doc: |-
  The identity operation. Useful for specifying data as a literal array, for scattering in a later step.

baseCommand: python3

requirements:
  InlineJavascriptRequirement: {}

inputs:
  input:
    label: The input data
    doc: |-
      The input data
    type: string
    format:
    - edam:format_2330

outputs:
  output:
    label: The output data
    doc: |-
      The output data
    type: string
    #format: edam:format_2330 # "'str' object does not support item assignment"
    outputBinding:
      outputEval: $(inputs.input)

#stdout: stdout

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
