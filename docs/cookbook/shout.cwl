cwlVersion: v1.2
class: CommandLineTool
doc: Write a copy of a text in capital letters.
requirements:
  DockerRequirement:
    dockerPull: debian:bookworm-slim
baseCommand: [sh, -c]
arguments:
- tr a-z A-Z < $(inputs.text.path) > shouted.txt
inputs:
  text:
    type: File
outputs:
  shouted:
    type: File
    outputBinding:
      glob: shouted.txt
