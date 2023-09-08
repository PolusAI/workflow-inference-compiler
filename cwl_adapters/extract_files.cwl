
#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: This class implements extracting array of files.

doc: |-
  This class implements extracting array of files

baseCommand: "true"
requirements:
  InlineJavascriptRequirement: {}

inputs:

  array:
    type:
      type: array
      items: File

  first_file:
    type: string?

  second_file:
    type: string?

outputs:

  first_file:
    label: output first file from array
    doc: output first file from array
    type: File
    outputBinding:
      outputEval: $(inputs.array[0])

  second_file:
    label: output second file from array
    doc: output second file from array
    type: File

    outputBinding:
      outputEval: $(inputs.array[1])

stdout: stdout

$namespaces:
  edam: https://edamontology.org/
