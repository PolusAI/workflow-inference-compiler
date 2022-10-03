#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: Downloads the NCI lgand database(s)

doc: |-
  Downloads the NCI lgand database(s)

baseCommand: cp
arguments: [$(inputs.database), $(runtime.outdir)]

requirements:
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerPull: jakefennick/scripts

inputs:
  database:
    type: string
    default: /NCIOpen.sdf # NOTE: Initial / required
    # /NCIOpen.sdf is mounted inside the docker image, . is outside of the image.

outputs:
  output_sdf_path:
    label: Path to the output sdf file
    doc: |-
      Path to the output sdf file
    type: File
    format: edam:format_3814 # sdf
    outputBinding:
      glob: $(inputs.database.slice(1)) # Remove initial /

stdout: stdout

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl