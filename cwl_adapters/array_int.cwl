#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: Generates an array of integers

doc: |
   Generates an array of integers

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

outputs:
  array:
    type: int[]
    outputBinding:
      outputEval: |
        ${
          var ints = [];
          for (var i = inputs.minval; i <= inputs.maxval; i += inputs.step) {
            ints.push(i);
          }
          return ints;
        }
