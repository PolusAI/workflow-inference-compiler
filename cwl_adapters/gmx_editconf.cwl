#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: Wrapper class for the GROMACS editconf module.

doc: |-
  The GROMACS solvate module generates a box around the selected structure.

baseCommand: gmx
arguments: ["-nobackup", "-nocopyright", "editconf"]

hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/biobb_md:3.7.2--pyhdfd78af_0

inputs:
  input_crd_path:
    label: Path to the input GRO file
    doc: |-
      Path to the input GRO file
      Type: string
      File type: input
      Accepted formats: gro, pdb
      Example file: https://github.com/bioexcel/biobb_md/raw/master/biobb_md/test/data/gromacs/editconf.gro
    type: File
    format:
    - edam:format_2033
    - edam:format_1476
    inputBinding:
      position: 1
      prefix: -f

  output_crd_path:
    label: Path to the output GRO file
    doc: |-
      Path to the output GRO file
      Type: string
      File type: output
      Accepted formats: pdb, gro
      Example file: https://github.com/bioexcel/biobb_md/raw/master/biobb_md/test/reference/gromacs/ref_editconf.gro
    type: string
    format:
    - edam:format_2033
    - edam:format_1476
    inputBinding:
      position: 2
      prefix: -o
#    default: system.pdb
#    default: system.gro
    default: system.g96

  distance_to_molecule:
    type: float
    inputBinding:
      position: 3
      prefix: -d
    default: 1.0

  box_type:
    type: string
    inputBinding:
      position: 4
      prefix: -bt
    default: cubic

  align_principal_axes:
    type: int? # Group number to align (0 == system)
    inputBinding:
      position: 5
      prefix: -princ

outputs:
  output_crd_path:
    label: Path to the output GRO file
    doc: |-
      Path to the output GRO file
    type: File
    outputBinding:
      glob: $(inputs.output_crd_path)
#    format: edam:format_1476
    format: edam:format_2033

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
