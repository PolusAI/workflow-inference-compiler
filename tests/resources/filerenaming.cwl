#!/usr/bin/env cwl-runner

cwlVersion: v1.1
class: CommandLineTool

requirements:
  DockerRequirement:
    dockerPull: polusai/file-renaming-plugin:0.2.0
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
  inpDir:
    type: Directory
    inputBinding:
      prefix: --inpDir
  mapDirectory:
    type: string?
    inputBinding:
      prefix: --mapDirectory
  outDir:
    type: Directory
    inputBinding:
      prefix: --outDir
  outFilePattern:
    type: string
    inputBinding:
      prefix: --outFilePattern

outputs:
  outDir:
    type: Directory
    outputBinding:
      glob: $(inputs.outDir.basename)
