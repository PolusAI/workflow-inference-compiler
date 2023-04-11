#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: Generates an array of integers, cast to strings.

doc: |
   Generates an array of integers, cast to strings.

baseCommand: 'true'

requirements:
  InlineJavascriptRequirement: {}

inputs: 
  minval:
    type: int
    default: 1

  maxval:
    type: int

  step:
    type: int
    default: 1

  ints_array:
    type: string?

outputs:
  ints_array:
    type: string[]
    outputBinding:
      outputEval: |
        ${
          var ints = [];
          for (var i = inputs.minval; i <= inputs.maxval; i += inputs.step) {
            ints.push(i.toString());
          }
          return ints;
        }
