#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: Generate a scatter plot

doc: |-
  Generate a scatter plot

baseCommand: python3
arguments: ["/scatter_plot.py"]

hints:
  DockerRequirement:
    dockerPull: jakefennick/scatter_plot

inputs:
#  script:
#    type: string
#    inputBinding:
#      position: 1
#    default: /scatter_plot.py

  xs:
    type:
      type: array
      items: float
    inputBinding:
      position: 2
      prefix: --xs

  ys:
    type:
      type: array
      items: float
    inputBinding:
      position: 3
      prefix: --ys

  ys2:
    type:
      type: array
      items: float
    inputBinding:
      position: 4
      prefix: --ys2

  output_png_path:
    label: Path to the output png file
    doc: |-
      Path to the output png file
    type: string
    format:
    - edam:format_3603 # png
    inputBinding:
      position: 1
      prefix: --output_png_path
    default: scatter.png

outputs:
  output_png_path:
    label: Path to the output png file
    doc: |-
      Path to the output png file
    type: File
    format: edam:format_3603 # png
    outputBinding:
      glob: $(inputs.output_png_path)

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl