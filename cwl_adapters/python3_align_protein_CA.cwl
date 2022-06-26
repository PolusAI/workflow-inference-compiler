#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: Align each frame of a trajectory w.r.t. the protein CA atoms in a reference structure.

doc: |-
   Align each frame of a trajectory w.r.t. the protein CA atoms in a reference structure.

baseCommand: python3

inputs:
  script:
    type: File
    format:
    - edam:format_2330 # 'Textual format'
    inputBinding:
      position: 1

  input_gro_path:
    label: Path to the reference coordinate file
    doc: |-
      Path to the reference coordinate file
    type: File
    format:
    - edam:format_2033 # Gromacs structure *.gro
    inputBinding:
      position: 2

  input_trr_path:
    label: Path to the input trajectory file
    doc: |-
      Path to the input trajectory file
    type: File
    format:
    - edam:format_3910 # Gromacs trajectory *.trr
    inputBinding:
      position: 3

  output_trr_path:
    label: Path to the output trajectory file
    doc: |-
      Path to the output trajectory file
    type: string
    format:
    - edam:format_3910 # Gromacs trajectory *.trr
    inputBinding:
      position: 4

outputs:
  output_trr_path:
    label: Path to the output trajectory file
    doc: |-
      Path to the output trajectory file
    type: File
    format: edam:format_3910 # Gromacs trajectory *.trr
    outputBinding:
      glob: $(inputs.output_trr_path)

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl