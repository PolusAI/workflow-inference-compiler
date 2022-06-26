#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: Wrapper of the GROMACS rms module for performing a Root Mean Square deviation
  (RMSd) analysis from a given GROMACS compatible trajectory.

doc: |-
  GROMACS rms compares two structures by computing the root mean square deviation (RMSD), the size-independent rho similarity parameter (rho) or the scaled rho (rhosc).

baseCommand: gmx
arguments: ["rms", "-fit", "none"] # rot+trans, translation, none

hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/biobb_analysis:3.7.0--pyhdfd78af_1

inputs:
  input_structure_path:
    label: Path to the input structure file
    doc: |-
      Path to the input structure file
      Type: string
      File type: input
      Accepted formats: tpr, gro, g96, pdb, brk, ent
      Example file: https://github.com/bioexcel/biobb_analysis/raw/master/biobb_analysis/test/data/gromacs/topology.tpr
    type: File
    format:
    - edam:format_2333
    - edam:format_2033
    - edam:format_2033
    - edam:format_1476
    - edam:format_2033
    - edam:format_1476
    inputBinding:
      position: 1
      prefix: -s

  input_traj_path:
    label: Path to the GROMACS trajectory file
    doc: |-
      Path to the GROMACS trajectory file
      Type: string
      File type: input
      Accepted formats: xtc, trr, cpt, gro, g96, pdb, tng
      Example file: https://github.com/bioexcel/biobb_analysis/raw/master/biobb_analysis/test/data/gromacs/trajectory.trr
    type: File
    format:
    - edam:format_3875
    - edam:format_3910
    - edam:format_2333
    - edam:format_2033
    - edam:format_2033
    - edam:format_1476
    - edam:format_3876
    inputBinding:
      position: 2
      prefix: -f

  output_xvg_path:
    label: Path to the XVG output file
    doc: |-
      Path to the XVG output file
      Type: string
      File type: output
      Accepted formats: xvg
      Example file: https://github.com/bioexcel/biobb_analysis/raw/master/biobb_analysis/test/reference/gromacs/ref_rms.xvg
    type: string
    format:
    - edam:format_2030
    inputBinding:
      position: 3
      prefix: -o
    default: system.xvg

  input_index_path:
    label: Path to the GROMACS index file
    doc: |-
      Path to the GROMACS index file
      Type: string
      File type: input
      Accepted formats: ndx
      Example file: https://github.com/bioexcel/biobb_analysis/raw/master/biobb_analysis/test/data/gromacs/index.ndx
    type: File?
    format:
    - edam:format_2033
    inputBinding:
      prefix: -n

outputs:
  output_xvg_path:
    label: Path to the XVG output file
    doc: |-
      Path to the XVG output file
    type: File
    outputBinding:
      glob: $(inputs.output_xvg_path)
    format: edam:format_2030

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
