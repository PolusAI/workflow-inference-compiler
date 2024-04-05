#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: Generates an array of floats

doc: |
   Generates an array of floats

baseCommand: 'true'

requirements:
  InlineJavascriptRequirement: {}

inputs:
  minval:
    type: float
    default: 1

  maxval:
    type: float

  step:
    type: float
    default: 1

  decimal_places:
    type: int
    default: 6
    label: Round to number of decimal places of precision.
    doc: Round to number of decimal places of precision.

outputs:
  array:
    type: float[]
    outputBinding:
      outputEval: |
        ${
          var arr = [];
          for (var num = inputs.minval; num <= inputs.maxval; num += inputs.step) {
            // Round to decimal_places
            arr.push(parseFloat(num.toFixed(inputs.decimal_places)));
          }
          return arr;
        }
