cwlVersion: v1.2
class: CommandLineTool
doc: Put a text and its word count together in one report.
requirements:
  DockerRequirement:
    dockerPull: debian:bookworm-slim
baseCommand: [sh, -c]
arguments:
- cat $(inputs.text.path) $(inputs.counts.path) > report.txt
inputs:
  text:
    type: File
  counts:
    type: File
outputs:
  report:
    type: File
    outputBinding:
      glob: report.txt
