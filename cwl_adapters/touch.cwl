cwlVersion: v1.0

class: CommandLineTool

requirements:
  DockerRequirement:
    dockerPull: docker.io/bash:4.4
  InlineJavascriptRequirement: {}

baseCommand: touch

inputs:
  filename:
    type: string
    inputBinding:
      position: 1

outputs:
  file:
    type: File
    outputBinding:
      glob: $(inputs.filename)

