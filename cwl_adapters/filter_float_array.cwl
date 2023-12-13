#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: Filter array of floats based on input boolean array

doc: |-
  Filter array of floats bases on input boolean array

baseCommand: 'true'

requirements:
  InlineJavascriptRequirement: {}

inputs:

  input_array:
    type: float[]

  input_bool_array:
    type: boolean[]

  output_array:
    type: string?

outputs:

  output_array:
    type: float[]
    # now we need to filter the input array based on the input bool array containing true, false
    outputBinding:
      glob:
      outputEval: $(inputs.input_array.filter((_, i) => inputs.input_bool_array[i]))

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl