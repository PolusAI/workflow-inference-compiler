
#!/usr/bin/env cwl-runner
cwlVersion: v1.2

class: CommandLineTool

label: return a subset from a array based on input indices

doc: |-
  This class return a subset from a array based on input indices
baseCommand: "true"

requirements:
  InlineJavascriptRequirement: {}

inputs:
  input_indices:
    label: input indices list
    type: int[]

  input_array:
    label: Path to the input array
    type: Any

outputs:
  output_array:
    label: Path to the output array (subset)
    doc: |-
      Path to the output array (subset)
    type: Any
    outputBinding:
      outputEval: |
        ${
          const data = [];
          for (const index of inputs.input_indices) {
            data.push(inputs.input_array[index]);
          }
          return data;
        }

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
