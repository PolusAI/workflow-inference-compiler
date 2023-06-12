#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

requirements:
  DockerRequirement:
    dockerPull: polusai/ome-converter-plugin:0.3.0
  InitialWorkDirRequirement:
    listing:
    - writable: true
      entry: $(inputs.outDir)
  InlineJavascriptRequirement: {}

inputs:
  fileExtension:
    type: string
    inputBinding:
      prefix: --fileExtension
  filePattern:
    type: string
    inputBinding:
      prefix: --filePattern
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
