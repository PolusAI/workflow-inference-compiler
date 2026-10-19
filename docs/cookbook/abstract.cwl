cwlVersion: v1.2
class: CommandLineTool
doc: Write an abstract from a summary; this lab's tools call a summary an abstract.
requirements:
  DockerRequirement:
    dockerPull: debian:bookworm-slim
baseCommand: cat
inputs:
  input_abstract:
    type: File
    inputBinding:
      position: 1
outputs:
  abstract:
    type: File
    outputBinding:
      glob: abstract.txt
stdout: abstract.txt
