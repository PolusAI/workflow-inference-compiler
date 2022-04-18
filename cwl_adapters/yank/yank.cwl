#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: Calculate binding free energy of receptor ligand systems using Yank.

doc: |-
  Calculate binding free energy of receptor ligand systems using Yank.

baseCommand: yank
arguments: [$(inputs.command)]

inputs:
  command:
    type: string
    format: edam:format_string

  yaml:
    label: Input YAML script
    type: File
    format: edam:format_3750
    inputBinding:
      prefix: -y

  input_receptor_path:
    label: Input receptor pdb file
    type: File
    format:
    - edam:format_1476

  input_ligand_path:
    label: Input ligand mol2 file
    type: File
    format:
    - edam:format_3816

outputs:
  output_log_path:
    label: Output log file
    type: File
    format: edam:format_2330 # 'Textual format'
    outputBinding:
      glob: output/experiments/experiments.log

  output_complex_netcdf_path:
    label: Binary output file for the complex phase
    type: File
    format: edam:format_3650
    outputBinding:
      glob: output/experiments/complex.nc

  output_complex_trailblaze_crd_path:
    label: Binary coordinates file for the complex phase
    type: File?
    format: edam:format_3878
    outputBinding:
      glob: output/experiments/trailblaze/complex/coordinates.dcd

  output_complex_trailblaze_protocol_path:
    label: YAML thermodynamic path protocol file for the complex phase
    type: File?
    format: edam:format_3750
    outputBinding:
      glob: output/experiments/trailblaze/complex/protocol.yaml

  output_complex_setup_crd_path:
    label: Output coordinates file for the complex phase (AMBER crd)
    type: File
    format: edam:format_3878
    outputBinding:
      glob: output/setup/systems/**/complex.inpcrd

  output_complex_setup_top_path:
    label: Output topology file for the complex phase (AMBER crd)
    type: File
    format: edam:format_3881
    outputBinding:
      glob: output/setup/systems/**/complex.prmtop

  output_solvent_netcdf_path:
    label: Binary output file for the solvent phase
    type: File
    format: edam:format_3650
    outputBinding:
      glob: output/experiments/solvent.nc

  output_solvent_trailblaze_crd_path:
    label: Binary coordinates file for the solvent phase
    type: File?
    format: edam:format_3878
    outputBinding:
      glob: output/experiments/trailblaze/solvent/coordinates.dcd

  output_solvent_trailblaze_protocol_path:
    label: YAML thermodynamic path protocol file for the solvent phase
    type: File?
    format: edam:format_3750
    outputBinding:
      glob: output/experiments/trailblaze/solvent/protocol.yaml

  output_solvent_setup_crd_path:
    label: Output coordinates file for the solvent phase (AMBER crd)
    type: File
    format: edam:format_3878
    outputBinding:
      glob: output/setup/systems/**/solvent.inpcrd

  output_solvent_setup_top_path:
    label: Output topology file for the solvent phase (AMBER crd)
    type: File
    format: edam:format_3881
    outputBinding:
      glob: output/setup/systems/**/solvent.prmtop

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl