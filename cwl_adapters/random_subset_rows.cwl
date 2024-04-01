
#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: return subset indices for array with sample selection using a random seed

doc: |-
  This class return subset indices for array with sample selection using a random seed

baseCommand: python3

hints:
  DockerRequirement:
    dockerPull: cyangnyu/random_subset_rows

requirements:
  InlineJavascriptRequirement: {}

inputs:
  script:
    type: string
    inputBinding:
      position: 1
    default: /random_subset_rows.py

  input_file:
    label: Path to the input file
    type: File
    inputBinding:
      position: 2
      prefix: --input_file

  num_of_samples:
    label: number of samples selected
    type: int
    inputBinding:
      position: 3
      prefix: --num_of_samples

  random_seed:
    label: random seed for the sample selected
    type: int
    inputBinding:
      position: 4
      prefix: --random_seed

  output_file:
    label: output file of subset indices
    type: string?
    inputBinding:
      position: 5
      prefix: --output_file
    default: "Index.txt"

outputs:
  output_file:
    label: output file of subset indices
    type: File
    outputBinding:
      glob: $(inputs.output_file)

  output_indices:
    label: indices list
    type: int[]
    outputBinding:
      glob: $(inputs.output_file)
      loadContents: true
      outputEval: |
        ${
          // The correct line should be of the form
          // 1
          // 11
          return self[0].contents.split("\n").map(str => parseInt(str));
        }

stdout: stdout

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
