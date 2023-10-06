#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: check_linear_fit

doc: |-
  check_linear_fit

baseCommand: python3
arguments: [$(inputs.script)]

requirements:
  DockerRequirement:
    dockerPull: cyangnyu/check_linear_fit

inputs:
  script:
    type: string
    default: /check_linear_fit.py

  xs:
    type:
      type: array
      items: float
    inputBinding:
      prefix: --xs

  ys:
    type:
      type: array
      items: float
    inputBinding:
      prefix: --ys

  tol_quad:
    type: float
    inputBinding:
      prefix: --tol_quad

  slope_min:
    type: float
    inputBinding:
      prefix: --slope_min

  slope_max:
    type: float
    inputBinding:
      prefix: --slope_max

  output:
    type: string?
outputs:
  output:
    type:
      type: array
      items: string
    outputBinding:
      glob: $(inputs.output)

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl