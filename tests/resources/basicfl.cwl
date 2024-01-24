#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
  DockerRequirement:
    dockerPull: polusai/basic-flatfield-estimation-plugin:2.1.0
  InitialWorkDirRequirement:
    listing:
    - writable: true
      entry: $(inputs.outDir)
  InlineJavascriptRequirement: {}

inputs:
  filePattern:
    type: string
    inputBinding:
      prefix: --filePattern
  getDarkfield:
    type: boolean
    inputBinding:
      prefix: --getDarkfield
  groupBy:
    type: string?
    inputBinding:
      prefix: --groupBy
  inpDir:
    type: Directory
    inputBinding:
      prefix: --inpDir
  outDir:
    type: Directory
    inputBinding:
      prefix: --outDir

outputs:
  outDir:
    type: Directory
    outputBinding:
      glob: $(inputs.outDir.basename)
