#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow
steps:
  step_1_gromacs_tutorial_setup:
    run: gromacs_tutorial_setup.cwl
    in:
      step_1_pdb_output_pdb_path_input: step_1_pdb_output_pdb_path_input
      step_1_pdb_config_input: step_1_pdb_config_input
      step_4_editconf_config_input: step_4_editconf_config_input
      step_6_grompp_config_input: step_6_grompp_config_input
      step_7_genion_config_input: step_7_genion_config_input
    out:
    - step_1_pdb_output_pdb_path
    - step_2_fix_side_chain_output_pdb_path
    - step_3_pdb2gmx_output_gro_path
    - step_3_pdb2gmx_output_top_zip_path
    - step_4_editconf_output_gro_path
    - step_5_solvate_output_gro_path
    - step_5_solvate_output_top_zip_path
    - step_6_grompp_output_tpr_path
    - step_7_genion_output_gro_path
    - step_7_genion_output_top_zip_path
  step_2_grompp:
    in:
      config: step_2_grompp_config_input
      output_tpr_path: step_2_grompp_output_tpr_path_input
      input_gro_path: step_1_gromacs_tutorial_setup/step_7_genion_output_gro_path
      input_top_zip_path: step_1_gromacs_tutorial_setup/step_7_genion_output_top_zip_path
    run: biobb/biobb_adapters/cwl/biobb_md/grompp.cwl
    out:
    - output_tpr_path
  step_3_mdrun:
    run: biobb/biobb_adapters/cwl/biobb_md/mdrun.cwl
    in:
      input_tpr_path: step_2_grompp/output_tpr_path
    out:
    - output_trr_path
    - output_gro_path
    - output_edr_path
    - output_log_path
    - output_xtc_path
    - output_cpt_path
    - output_dhdl_path
  step_4_gmx_energy:
    in:
      config: step_4_gmx_energy_config_input
      output_xvg_path: step_4_gmx_energy_output_xvg_path_input
      input_energy_path: step_3_mdrun/output_edr_path
    run: biobb/biobb_adapters/cwl/biobb_analysis/gmx_energy.cwl
    out:
    - output_xvg_path
  step_5_grompp:
    in:
      config: step_5_grompp_config_input
      input_gro_path: step_3_mdrun/output_gro_path
      input_top_zip_path: step_1_gromacs_tutorial_setup/step_7_genion_output_top_zip_path
    run: biobb/biobb_adapters/cwl/biobb_md/grompp.cwl
    out:
    - output_tpr_path
  step_6_mdrun:
    run: biobb/biobb_adapters/cwl/biobb_md/mdrun.cwl
    in:
      input_tpr_path: step_5_grompp/output_tpr_path
    out:
    - output_trr_path
    - output_gro_path
    - output_edr_path
    - output_log_path
    - output_xtc_path
    - output_cpt_path
    - output_dhdl_path
  step_7_gmx_energy:
    in:
      config: step_7_gmx_energy_config_input
      output_xvg_path: step_7_gmx_energy_output_xvg_path_input
      input_energy_path: step_6_mdrun/output_edr_path
    run: biobb/biobb_adapters/cwl/biobb_analysis/gmx_energy.cwl
    out:
    - output_xvg_path
  step_8_grompp:
    in:
      config: step_8_grompp_config_input
      input_gro_path: step_6_mdrun/output_gro_path
      input_top_zip_path: step_1_gromacs_tutorial_setup/step_7_genion_output_top_zip_path
    run: biobb/biobb_adapters/cwl/biobb_md/grompp.cwl
    out:
    - output_tpr_path
  step_9_mdrun:
    run: biobb/biobb_adapters/cwl/biobb_md/mdrun.cwl
    in:
      input_tpr_path: step_8_grompp/output_tpr_path
    out:
    - output_trr_path
    - output_gro_path
    - output_edr_path
    - output_log_path
    - output_xtc_path
    - output_cpt_path
    - output_dhdl_path
  step_10_gmx_energy:
    in:
      config: step_10_gmx_energy_config_input
      output_xvg_path: step_10_gmx_energy_output_xvg_path_input
      input_energy_path: step_9_mdrun/output_edr_path
    run: biobb/biobb_adapters/cwl/biobb_analysis/gmx_energy.cwl
    out:
    - output_xvg_path
  step_11_grompp:
    in:
      config: step_11_grompp_config_input
      output_tpr_path: step_11_grompp_output_tpr_path_input
      input_gro_path: step_9_mdrun/output_gro_path
      input_top_zip_path: step_1_gromacs_tutorial_setup/step_7_genion_output_top_zip_path
    run: biobb/biobb_adapters/cwl/biobb_md/grompp.cwl
    out:
    - output_tpr_path
  step_12_mdrun:
    in:
      output_gro_path: step_12_mdrun_output_gro_path_input
      input_tpr_path: step_11_grompp/output_tpr_path
    run: biobb/biobb_adapters/cwl/biobb_md/mdrun.cwl
    out:
    - output_trr_path
    - output_gro_path
    - output_edr_path
    - output_log_path
    - output_xtc_path
    - output_cpt_path
    - output_dhdl_path
  step_13_gmx_rms:
    in:
      config: step_13_gmx_rms_config_input
      output_xvg_path: step_13_gmx_rms_output_xvg_path_input
      input_structure_path: step_11_grompp/output_tpr_path
      input_traj_path: step_12_mdrun/output_trr_path
    run: biobb/biobb_adapters/cwl/biobb_analysis/gmx_rms.cwl
    out:
    - output_xvg_path
  step_14_gmx_rms:
    in:
      config: step_14_gmx_rms_config_input
      output_xvg_path: step_14_gmx_rms_output_xvg_path_input
      input_structure_path: step_2_grompp/output_tpr_path
      input_traj_path: step_12_mdrun/output_trr_path
    run: biobb/biobb_adapters/cwl/biobb_analysis/gmx_rms.cwl
    out:
    - output_xvg_path
  step_15_gmx_rgyr:
    in:
      config: step_15_gmx_rgyr_config_input
      input_structure_path: step_2_grompp/output_tpr_path
      input_traj_path: step_12_mdrun/output_trr_path
    run: biobb/biobb_adapters/cwl/biobb_analysis/gmx_rgyr.cwl
    out:
    - output_xvg_path
  step_16_gmx_image:
    in:
      config: step_16_gmx_image_config_input
      input_top_path: step_11_grompp/output_tpr_path
      input_traj_path: step_12_mdrun/output_trr_path
    run: biobb/biobb_adapters/cwl/biobb_analysis/gmx_image.cwl
    out:
    - output_traj_path
  step_17_gmx_trjconv_str:
    in:
      config: step_17_gmx_trjconv_str_config_input
      input_top_path: step_11_grompp/output_tpr_path
      input_structure_path: step_12_mdrun/output_gro_path
    run: biobb/biobb_adapters/cwl/biobb_analysis/gmx_trjconv_str.cwl
    out:
    - output_str_path
requirements:
  SubworkflowFeatureRequirement:
    class: SubworkflowFeatureRequirement
inputs:
  step_1_pdb_output_pdb_path_input: string
  step_1_pdb_config_input: string
  step_4_editconf_config_input: string
  step_6_grompp_config_input: string
  step_7_genion_config_input: string
  step_2_grompp_config_input: string
  step_2_grompp_output_tpr_path_input: string
  step_4_gmx_energy_config_input: string
  step_4_gmx_energy_output_xvg_path_input: string
  step_5_grompp_config_input: string
  step_7_gmx_energy_config_input: string
  step_7_gmx_energy_output_xvg_path_input: string
  step_8_grompp_config_input: string
  step_10_gmx_energy_config_input: string
  step_10_gmx_energy_output_xvg_path_input: string
  step_11_grompp_config_input: string
  step_11_grompp_output_tpr_path_input: string
  step_12_mdrun_output_gro_path_input: string
  step_13_gmx_rms_config_input: string
  step_13_gmx_rms_output_xvg_path_input: string
  step_14_gmx_rms_config_input: string
  step_14_gmx_rms_output_xvg_path_input: string
  step_15_gmx_rgyr_config_input: string
  step_16_gmx_image_config_input: string
  step_17_gmx_trjconv_str_config_input: string
outputs:
  step_1_gromacs_tutorial_setup_step_1_pdb_output_pdb_path:
    type: File
    outputSource: step_1_gromacs_tutorial_setup/step_1_pdb_output_pdb_path
  step_1_gromacs_tutorial_setup_step_2_fix_side_chain_output_pdb_path:
    type: File
    outputSource: step_1_gromacs_tutorial_setup/step_2_fix_side_chain_output_pdb_path
  step_1_gromacs_tutorial_setup_step_3_pdb2gmx_output_gro_path:
    type: File
    outputSource: step_1_gromacs_tutorial_setup/step_3_pdb2gmx_output_gro_path
  step_1_gromacs_tutorial_setup_step_3_pdb2gmx_output_top_zip_path:
    type: File
    outputSource: step_1_gromacs_tutorial_setup/step_3_pdb2gmx_output_top_zip_path
  step_1_gromacs_tutorial_setup_step_4_editconf_output_gro_path:
    type: File
    outputSource: step_1_gromacs_tutorial_setup/step_4_editconf_output_gro_path
  step_1_gromacs_tutorial_setup_step_5_solvate_output_gro_path:
    type: File
    outputSource: step_1_gromacs_tutorial_setup/step_5_solvate_output_gro_path
  step_1_gromacs_tutorial_setup_step_5_solvate_output_top_zip_path:
    type: File
    outputSource: step_1_gromacs_tutorial_setup/step_5_solvate_output_top_zip_path
  step_1_gromacs_tutorial_setup_step_6_grompp_output_tpr_path:
    type: File
    outputSource: step_1_gromacs_tutorial_setup/step_6_grompp_output_tpr_path
  step_1_gromacs_tutorial_setup_step_7_genion_output_gro_path:
    type: File
    outputSource: step_1_gromacs_tutorial_setup/step_7_genion_output_gro_path
  step_1_gromacs_tutorial_setup_step_7_genion_output_top_zip_path:
    type: File
    outputSource: step_1_gromacs_tutorial_setup/step_7_genion_output_top_zip_path
  step_2_grompp_output_tpr_path:
    type: File
    outputSource: step_2_grompp/output_tpr_path
  step_3_mdrun_output_trr_path:
    type: File
    outputSource: step_3_mdrun/output_trr_path
  step_3_mdrun_output_gro_path:
    type: File
    outputSource: step_3_mdrun/output_gro_path
  step_3_mdrun_output_edr_path:
    type: File
    outputSource: step_3_mdrun/output_edr_path
  step_3_mdrun_output_log_path:
    type: File
    outputSource: step_3_mdrun/output_log_path
  step_3_mdrun_output_cpt_path:
    type: File
    outputSource: step_3_mdrun/output_cpt_path
  step_4_gmx_energy_output_xvg_path:
    type: File
    outputSource: step_4_gmx_energy/output_xvg_path
  step_5_grompp_output_tpr_path:
    type: File
    outputSource: step_5_grompp/output_tpr_path
  step_6_mdrun_output_trr_path:
    type: File
    outputSource: step_6_mdrun/output_trr_path
  step_6_mdrun_output_gro_path:
    type: File
    outputSource: step_6_mdrun/output_gro_path
  step_6_mdrun_output_edr_path:
    type: File
    outputSource: step_6_mdrun/output_edr_path
  step_6_mdrun_output_log_path:
    type: File
    outputSource: step_6_mdrun/output_log_path
  step_6_mdrun_output_cpt_path:
    type: File
    outputSource: step_6_mdrun/output_cpt_path
  step_7_gmx_energy_output_xvg_path:
    type: File
    outputSource: step_7_gmx_energy/output_xvg_path
  step_8_grompp_output_tpr_path:
    type: File
    outputSource: step_8_grompp/output_tpr_path
  step_9_mdrun_output_trr_path:
    type: File
    outputSource: step_9_mdrun/output_trr_path
  step_9_mdrun_output_gro_path:
    type: File
    outputSource: step_9_mdrun/output_gro_path
  step_9_mdrun_output_edr_path:
    type: File
    outputSource: step_9_mdrun/output_edr_path
  step_9_mdrun_output_log_path:
    type: File
    outputSource: step_9_mdrun/output_log_path
  step_9_mdrun_output_cpt_path:
    type: File
    outputSource: step_9_mdrun/output_cpt_path
  step_10_gmx_energy_output_xvg_path:
    type: File
    outputSource: step_10_gmx_energy/output_xvg_path
  step_11_grompp_output_tpr_path:
    type: File
    outputSource: step_11_grompp/output_tpr_path
  step_12_mdrun_output_trr_path:
    type: File
    outputSource: step_12_mdrun/output_trr_path
  step_12_mdrun_output_gro_path:
    type: File
    outputSource: step_12_mdrun/output_gro_path
  step_12_mdrun_output_edr_path:
    type: File
    outputSource: step_12_mdrun/output_edr_path
  step_12_mdrun_output_log_path:
    type: File
    outputSource: step_12_mdrun/output_log_path
  step_12_mdrun_output_cpt_path:
    type: File
    outputSource: step_12_mdrun/output_cpt_path
  step_13_gmx_rms_output_xvg_path:
    type: File
    outputSource: step_13_gmx_rms/output_xvg_path
  step_14_gmx_rms_output_xvg_path:
    type: File
    outputSource: step_14_gmx_rms/output_xvg_path
  step_15_gmx_rgyr_output_xvg_path:
    type: File
    outputSource: step_15_gmx_rgyr/output_xvg_path
  step_16_gmx_image_output_traj_path:
    type: File
    outputSource: step_16_gmx_image/output_traj_path
  step_17_gmx_trjconv_str_output_str_path:
    type: File
    outputSource: step_17_gmx_trjconv_str/output_str_path
