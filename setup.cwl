#!/usr/bin/env cwl-runner
steps:
  setup_step_1_pdb:
    in:
      output_pdb_path: setup_step_1_pdb___output_pdb_path
      config: setup_step_1_pdb___config
    run: biobb/biobb_adapters/cwl/biobb_io/pdb.cwl
    out:
    - output_pdb_path
  setup_step_2_fix_side_chain:
    run: biobb/biobb_adapters/cwl/biobb_model/fix_side_chain.cwl
    in:
      input_pdb_path: setup_step_1_pdb/output_pdb_path
    out:
    - output_pdb_path
  setup_step_3_pdb2gmx:
    run: biobb/biobb_adapters/cwl/biobb_md/pdb2gmx.cwl
    in:
      input_pdb_path: setup_step_2_fix_side_chain/output_pdb_path
    out:
    - output_gro_path
    - output_top_zip_path
  setup_step_4_editconf:
    in:
      config: setup_step_4_editconf___config
      input_gro_path: setup_step_3_pdb2gmx/output_gro_path
    run: biobb/biobb_adapters/cwl/biobb_md/editconf.cwl
    out:
    - output_gro_path
  setup_step_5_solvate:
    run: biobb/biobb_adapters/cwl/biobb_md/solvate.cwl
    in:
      input_solute_gro_path: setup_step_4_editconf/output_gro_path
      input_top_zip_path: setup_step_3_pdb2gmx/output_top_zip_path
    out:
    - output_gro_path
    - output_top_zip_path
  setup_step_6_grompp:
    in:
      config: setup_step_6_grompp___config
      input_gro_path: setup_step_5_solvate/output_gro_path
      input_top_zip_path: setup_step_5_solvate/output_top_zip_path
    run: biobb/biobb_adapters/cwl/biobb_md/grompp.cwl
    out:
    - output_tpr_path
  setup_step_7_genion:
    in:
      config: setup_step_7_genion___config
      input_tpr_path: setup_step_6_grompp/output_tpr_path
      input_top_zip_path: setup_step_5_solvate/output_top_zip_path
    run: biobb/biobb_adapters/cwl/biobb_md/genion.cwl
    out:
    - output_gro_path
    - output_top_zip_path
cwlVersion: v1.0
class: Workflow
inputs:
  setup_step_1_pdb___output_pdb_path: string
  setup_step_1_pdb___config: string
  setup_step_4_editconf___config: string
  setup_step_6_grompp___config: string
  setup_step_7_genion___config: string
outputs:
  setup_step_1_pdb___output_pdb_path:
    type: File
    outputSource: setup_step_1_pdb/output_pdb_path
  setup_step_2_fix_side_chain___output_pdb_path:
    type: File
    outputSource: setup_step_2_fix_side_chain/output_pdb_path
  setup_step_3_pdb2gmx___output_gro_path:
    type: File
    outputSource: setup_step_3_pdb2gmx/output_gro_path
  setup_step_3_pdb2gmx___output_top_zip_path:
    type: File
    outputSource: setup_step_3_pdb2gmx/output_top_zip_path
  setup_step_4_editconf___output_gro_path:
    type: File
    outputSource: setup_step_4_editconf/output_gro_path
  setup_step_5_solvate___output_gro_path:
    type: File
    outputSource: setup_step_5_solvate/output_gro_path
  setup_step_5_solvate___output_top_zip_path:
    type: File
    outputSource: setup_step_5_solvate/output_top_zip_path
  setup_step_6_grompp___output_tpr_path:
    type: File
    outputSource: setup_step_6_grompp/output_tpr_path
  setup_step_7_genion___output_gro_path:
    type: File
    outputSource: setup_step_7_genion/output_gro_path
  setup_step_7_genion___output_top_zip_path:
    type: File
    outputSource: setup_step_7_genion/output_top_zip_path
