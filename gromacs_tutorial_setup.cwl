#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow
steps:
  step_1_pdb:
    in:
      output_pdb_path: step_1_pdb_output_pdb_path_input
      config: step_1_pdb_config_input
    run: biobb/biobb_adapters/cwl/biobb_io/pdb.cwl
    out:
    - output_pdb_path
  step_2_fix_side_chain:
    run: biobb/biobb_adapters/cwl/biobb_model/fix_side_chain.cwl
    in:
      input_pdb_path: step_1_pdb/output_pdb_path
    out:
    - output_pdb_path
  step_3_pdb2gmx:
    run: biobb/biobb_adapters/cwl/biobb_md/pdb2gmx.cwl
    in:
      input_pdb_path: step_2_fix_side_chain/output_pdb_path
    out:
    - output_gro_path
    - output_top_zip_path
  step_4_editconf:
    in:
      config: step_4_editconf_config_input
      input_gro_path: step_3_pdb2gmx/output_gro_path
    run: biobb/biobb_adapters/cwl/biobb_md/editconf.cwl
    out:
    - output_gro_path
  step_5_solvate:
    run: biobb/biobb_adapters/cwl/biobb_md/solvate.cwl
    in:
      input_solute_gro_path: step_4_editconf/output_gro_path
      input_top_zip_path: step_3_pdb2gmx/output_top_zip_path
    out:
    - output_gro_path
    - output_top_zip_path
  step_6_grompp:
    in:
      config: step_6_grompp_config_input
      input_gro_path: step_5_solvate/output_gro_path
      input_top_zip_path: step_5_solvate/output_top_zip_path
    run: biobb/biobb_adapters/cwl/biobb_md/grompp.cwl
    out:
    - output_tpr_path
  step_7_genion:
    in:
      config: step_7_genion_config_input
      input_tpr_path: step_6_grompp/output_tpr_path
      input_top_zip_path: step_5_solvate/output_top_zip_path
    run: biobb/biobb_adapters/cwl/biobb_md/genion.cwl
    out:
    - output_gro_path
    - output_top_zip_path
inputs:
  step_1_pdb_output_pdb_path_input: string
  step_1_pdb_config_input: string
  step_4_editconf_config_input: string
  step_6_grompp_config_input: string
  step_7_genion_config_input: string
outputs:
  step_1_pdb_output_pdb_path:
    type: File
    outputSource: step_1_pdb/output_pdb_path
  step_2_fix_side_chain_output_pdb_path:
    type: File
    outputSource: step_2_fix_side_chain/output_pdb_path
  step_3_pdb2gmx_output_gro_path:
    type: File
    outputSource: step_3_pdb2gmx/output_gro_path
  step_3_pdb2gmx_output_top_zip_path:
    type: File
    outputSource: step_3_pdb2gmx/output_top_zip_path
  step_4_editconf_output_gro_path:
    type: File
    outputSource: step_4_editconf/output_gro_path
  step_5_solvate_output_gro_path:
    type: File
    outputSource: step_5_solvate/output_gro_path
  step_5_solvate_output_top_zip_path:
    type: File
    outputSource: step_5_solvate/output_top_zip_path
  step_6_grompp_output_tpr_path:
    type: File
    outputSource: step_6_grompp/output_tpr_path
  step_7_genion_output_gro_path:
    type: File
    outputSource: step_7_genion/output_gro_path
  step_7_genion_output_top_zip_path:
    type: File
    outputSource: step_7_genion/output_top_zip_path
