cwlVersion: v1.2
class: CommandLineTool
doc: Write a report from a summary.
requirements:
  DockerRequirement:
    dockerPull: debian:bookworm-slim
baseCommand: cat
inputs:
  input_summary:
    type: File
    inputBinding:
      position: 1
outputs:
  report:
    type: File
    outputBinding:
      glob: report.txt
stdout: report.txt
