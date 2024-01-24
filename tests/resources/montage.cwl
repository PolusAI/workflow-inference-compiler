class: CommandLineTool
cwlVersion: v1.2
inputs:
  filePattern:
    inputBinding:
      prefix: --filePattern
    type: string
  flipAxis:
    inputBinding:
      prefix: --flipAxis
    type: string?
  gridSpacing:
    inputBinding:
      prefix: --gridSpacing
    type: string?
  imageSpacing:
    inputBinding:
      prefix: --imageSpacing
    type: string?
  inpDir:
    inputBinding:
      prefix: --inpDir
    type: Directory
  layout:
    inputBinding:
      prefix: --layout
    type: string?
  outDir:
    inputBinding:
      prefix: --outDir
    type: Directory
outputs:
  outDir:
    outputBinding:
      glob: $(inputs.outDir.basename)
    type: Directory
requirements:
  DockerRequirement:
    dockerPull: polusai/montage-plugin:0.5.0
  InitialWorkDirRequirement:
    listing:
    - entry: $(inputs.outDir)
      writable: true
  InlineJavascriptRequirement: {}
