#!/usr/bin/env cwl-runner
cwlVersion: v1.2

class: CommandLineTool

label: Filter array based on input boolean array

doc: |-
  Filter array of Files based on input boolean array

baseCommand: 'true'

requirements:
  InlineJavascriptRequirement: {}

inputs:

  input_array:
    type: Any[]

  input_bool_array:
    type: boolean[]

outputs:

  output_array:
    type: Any[]
    # now we need to filter the input array based on the input bool array containing true, false
    outputBinding:
      glob:
      outputEval: $(inputs.input_array.filter((_, i) => inputs.input_bool_array[i]))

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl