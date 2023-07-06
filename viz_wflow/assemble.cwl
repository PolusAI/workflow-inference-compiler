#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

requirements:
  DockerRequirement:
    dockerPull: polusai/image-assembler-plugin:1.3.0-dev0
  InitialWorkDirRequirement:
    listing:
    - writable: true
      entry: $(inputs.outDir)
  InlineJavascriptRequirement: {}

inputs:
  imgPath:
    type: Directory
    inputBinding:
      prefix: --imgPath
  outDir:
    type: Directory
    inputBinding:
      prefix: --outDir
  preview:
    type: boolean?
    inputBinding:
      prefix: --preview
  stitchPath:
    type: Directory
    inputBinding:
      prefix: --stitchPath
  timesliceNaming:
    type: boolean?
    inputBinding:
      prefix: --timesliceNaming

outputs:
  outDir:
    type: Directory
    outputBinding:
      glob: $(inputs.outDir.basename)
