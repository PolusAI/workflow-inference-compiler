
#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: This class implements zipping two 1D arrays of Files together

doc: |-
  This class implements zipping two 1D arrays of Files together

baseCommand: "true"
requirements:
  - class: InlineJavascriptRequirement

inputs:

  first_array:
    type:
      type: array
      items: File

  second_array:
    type:
      type: array
      items: File

  zipped_array:
    type: string?

outputs:

  zipped_array:
    label: output array of array
    doc: output array of array
    type:
      type: array
      items:
        type: array
        items: File

    outputBinding:
      outputEval: |
        ${
          var lst = [];
          for(var i = 0; i < inputs.first_array.length; i++)
          {
              lst.push([inputs.first_array[i], inputs.second_array[i]]);
          }
          return lst;
        }

stdout: stdout

$namespaces:
  edam: https://edamontology.org/