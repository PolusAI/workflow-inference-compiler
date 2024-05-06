
#!/usr/bin/env cwl-runner
cwlVersion: v1.2

class: CommandLineTool

label: Duplicates an input n times

doc: |-
  Duplicates an input n times

baseCommand: "true"

requirements:
  InlineJavascriptRequirement: {}

inputs:
  input_scalar:
    type: Any

  output_array_length:
    type: int

outputs:
  output_array:
    type: Any[]
    outputBinding:
      outputEval: |
        ${
          var any_array = [];
          for (var i = 0; i < inputs.output_array_length; i++) {
            any_array.push(inputs.input_scalar)
          }
          return any_array;
        }

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
