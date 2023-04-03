cwlVersion: v1.0

class: CommandLineTool

baseCommand: echo

inputs:
  message:
    type: string
    #default: "Hello World"
    inputBinding:
      position: 1

outputs:
  stdout:
    type: File
    outputBinding:
      glob: stdout

stdout: stdout