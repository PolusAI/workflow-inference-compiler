cwlVersion: v1.0

class: CommandLineTool

baseCommand: echo
arguments: ["$(inputs.message1) $(inputs.message2) $(inputs.message3)"]

inputs:
  message1:
    type: string

  message2:
    type: string

  message3:
    type: string

outputs:
  stdout:
    type: File
    outputBinding:
      glob: stdout

stdout: stdout