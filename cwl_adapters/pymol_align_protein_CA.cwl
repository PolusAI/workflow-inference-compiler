#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: Run a pymol script

doc: |-
  Run a pymol script

baseCommand: pymol
arguments: ["-rcQ", $(inputs.script), "--", $(inputs.input_1_path), $(inputs.input_2_path), $(inputs.input_3_path), $(inputs.input_4_path), $(inputs.output_file_path)]
# NOTE: Based on the last example here
# See https://pymolwiki.org/index.php/Command_Line_Options

hints:
  DockerRequirement:
    dockerPull: jakefennick/scripts

inputs:
  script:
    type: string
    inputBinding:
      position: 1
    default: /align_protein_CA_pymol.py # NOTE: Initial / required

  input_1_path:
    type: File
    format:
    - edam:format_1476 # pdb

  input_2_path:
    type: File
    format:
    - edam:format_1476 # pdb

  input_3_path:
    type: File
    format:
    - edam:format_2033 # Gromacs structure *.gro

  input_4_path:
    type: File
    format:
    - edam:format_3910 # Gromacs trajectory *.trr

  output_file_path:
    label: Path to the output file
    doc: |-
      Path to the output file
    type: string
    format:
    - edam:format_1476 # pdb
    default: system.pdb

outputs:
  output_file_path:
    label: Path to the output file
    doc: |-
      Path to the output file
    type: File
    format: edam:format_1476 # pdb
    outputBinding:
      glob: $(inputs.output_file_path)

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl