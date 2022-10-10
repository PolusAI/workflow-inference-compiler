#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: Downloads the SMACC ligand database(s) https://smacc.mml.unc.edu

doc: |-
  Downloads the SMACC ligand database(s) https://smacc.mml.unc.edu

baseCommand: cp
arguments: [$(inputs.database), $(runtime.outdir)]

requirements:
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerPull: jakefennick/data

inputs:
  database:
    type: string
    default: /ncats_target_based_curated.xlsx # NOTE: Initial / required
    # /ncats_target_based_curated.xlsx is mounted inside the docker image, . is outside of the image.

outputs:
  output_excel_path:
    label: Path to the output xlsx file
    doc: |-
      Path to the output xlsx file
    type: File
    format: edam:format_3620
    outputBinding:
      glob: $(inputs.database.slice(1)) # Remove initial /

stdout: stdout

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl