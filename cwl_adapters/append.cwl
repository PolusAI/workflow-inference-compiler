class: CommandLineTool
cwlVersion: v1.0

requirements:
  ShellCommandRequirement: {}
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
    - $(inputs.file)

# NOTE: This mutates file in-place!

inputs:
  str:
    type: string
    inputBinding:
      shellQuote: false
      position: 1
      prefix: echo

  file:
    type: File
    inputBinding:
      shellQuote: false
      position: 2
      prefix: ">>"
      #prefix: "| tee --append"

outputs:
  file:
    type: File
    outputBinding:
      glob: $(inputs.file.basename)
