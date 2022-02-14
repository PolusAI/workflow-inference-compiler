#!/usr/bin/env cwl-runner
steps:
  wf_step_1_setup.yml:
    run: setup.cwl
    in:
      setup_step_1_pdb___output_pdb_path: wf_step_1_setup.yml___setup_step_1_pdb___output_pdb_path
      setup_step_1_pdb___config: wf_step_1_setup.yml___setup_step_1_pdb___config
      setup_step_4_editconf___config: wf_step_1_setup.yml___setup_step_4_editconf___config
      setup_step_6_grompp___config: wf_step_1_setup.yml___setup_step_6_grompp___config
      setup_step_7_genion___config: wf_step_1_setup.yml___setup_step_7_genion___config
    out:
    - setup_step_1_pdb___output_pdb_path
    - setup_step_2_fix_side_chain___output_pdb_path
    - setup_step_3_pdb2gmx___output_gro_path
    - setup_step_3_pdb2gmx___output_top_zip_path
    - setup_step_4_editconf___output_gro_path
    - setup_step_5_solvate___output_gro_path
    - setup_step_5_solvate___output_top_zip_path
    - setup_step_6_grompp___output_tpr_path
    - setup_step_7_genion___output_gro_path
    - setup_step_7_genion___output_top_zip_path
  wf_step_2_min.yml:
    run: min.cwl
    in:
      min_step_1_grompp___config: wf_step_2_min.yml___min_step_1_grompp___config
      min_step_1_grompp___output_tpr_path: wf_step_2_min.yml___min_step_1_grompp___output_tpr_path
      min_step_1_grompp___input_gro_path: wf_step_1_setup.yml/setup_step_7_genion___output_gro_path
      min_step_1_grompp___input_top_zip_path: wf_step_1_setup.yml/setup_step_7_genion___output_top_zip_path
      min_step_3_gmx_energy___config: wf_step_2_min.yml___min_step_3_gmx_energy___config
      min_step_3_gmx_energy___output_xvg_path: wf_step_2_min.yml___min_step_3_gmx_energy___output_xvg_path
    out:
    - min_step_1_grompp___output_tpr_path
    - min_step_2_mdrun___output_trr_path
    - min_step_2_mdrun___output_gro_path
    - min_step_2_mdrun___output_edr_path
    - min_step_2_mdrun___output_log_path
    - min_step_2_mdrun___output_cpt_path
    - min_step_3_gmx_energy___output_xvg_path
  wf_step_3_nvt.yml:
    run: nvt.cwl
    in:
      nvt_step_1_grompp___config: wf_step_3_nvt.yml___nvt_step_1_grompp___config
      nvt_step_1_grompp___input_gro_path: wf_step_2_min.yml/min_step_2_mdrun___output_gro_path
      nvt_step_1_grompp___input_top_zip_path: wf_step_1_setup.yml/setup_step_7_genion___output_top_zip_path
      nvt_step_3_gmx_energy___config: wf_step_3_nvt.yml___nvt_step_3_gmx_energy___config
      nvt_step_3_gmx_energy___output_xvg_path: wf_step_3_nvt.yml___nvt_step_3_gmx_energy___output_xvg_path
    out:
    - nvt_step_1_grompp___output_tpr_path
    - nvt_step_2_mdrun___output_trr_path
    - nvt_step_2_mdrun___output_gro_path
    - nvt_step_2_mdrun___output_edr_path
    - nvt_step_2_mdrun___output_log_path
    - nvt_step_2_mdrun___output_cpt_path
    - nvt_step_3_gmx_energy___output_xvg_path
  wf_step_4_npt.yml:
    run: npt.cwl
    in:
      npt_step_1_grompp___config: wf_step_4_npt.yml___npt_step_1_grompp___config
      npt_step_1_grompp___input_gro_path: wf_step_3_nvt.yml/nvt_step_2_mdrun___output_gro_path
      npt_step_1_grompp___input_top_zip_path: wf_step_1_setup.yml/setup_step_7_genion___output_top_zip_path
      npt_step_3_gmx_energy___config: wf_step_4_npt.yml___npt_step_3_gmx_energy___config
      npt_step_3_gmx_energy___output_xvg_path: wf_step_4_npt.yml___npt_step_3_gmx_energy___output_xvg_path
    out:
    - npt_step_1_grompp___output_tpr_path
    - npt_step_2_mdrun___output_trr_path
    - npt_step_2_mdrun___output_gro_path
    - npt_step_2_mdrun___output_edr_path
    - npt_step_2_mdrun___output_log_path
    - npt_step_2_mdrun___output_cpt_path
    - npt_step_3_gmx_energy___output_xvg_path
  wf_step_5_grompp:
    in:
      config: wf_step_5_grompp___config
      output_tpr_path: wf_step_5_grompp___output_tpr_path
      input_gro_path: wf_step_4_npt.yml/npt_step_2_mdrun___output_gro_path
      input_top_zip_path: wf_step_1_setup.yml/setup_step_7_genion___output_top_zip_path
    run: biobb/biobb_adapters/cwl/biobb_md/grompp.cwl
    out:
    - output_tpr_path
  wf_step_6_mdrun:
    in:
      output_gro_path: wf_step_6_mdrun___output_gro_path
      input_tpr_path: wf_step_5_grompp/output_tpr_path
    run: biobb/biobb_adapters/cwl/biobb_md/mdrun.cwl
    out:
    - output_trr_path
    - output_gro_path
    - output_edr_path
    - output_log_path
    - output_xtc_path
    - output_cpt_path
    - output_dhdl_path
  wf_step_7_gmx_rms:
    in:
      config: wf_step_7_gmx_rms___config
      output_xvg_path: wf_step_7_gmx_rms___output_xvg_path
      input_structure_path: wf_step_5_grompp/output_tpr_path
      input_traj_path: wf_step_6_mdrun/output_trr_path
    run: biobb/biobb_adapters/cwl/biobb_analysis/gmx_rms.cwl
    out:
    - output_xvg_path
  wf_step_8_gmx_rms:
    in:
      config: wf_step_8_gmx_rms___config
      output_xvg_path: wf_step_8_gmx_rms___output_xvg_path
      input_structure_path: wf_step_2_min.yml/min_step_1_grompp___output_tpr_path
      input_traj_path: wf_step_6_mdrun/output_trr_path
    run: biobb/biobb_adapters/cwl/biobb_analysis/gmx_rms.cwl
    out:
    - output_xvg_path
  wf_step_9_gmx_rgyr:
    in:
      config: wf_step_9_gmx_rgyr___config
      input_structure_path: wf_step_2_min.yml/min_step_1_grompp___output_tpr_path
      input_traj_path: wf_step_6_mdrun/output_trr_path
    run: biobb/biobb_adapters/cwl/biobb_analysis/gmx_rgyr.cwl
    out:
    - output_xvg_path
  wf_step_10_gmx_image:
    in:
      config: wf_step_10_gmx_image___config
      input_top_path: wf_step_5_grompp/output_tpr_path
      input_traj_path: wf_step_6_mdrun/output_trr_path
    run: biobb/biobb_adapters/cwl/biobb_analysis/gmx_image.cwl
    out:
    - output_traj_path
  wf_step_11_gmx_trjconv_str:
    in:
      config: wf_step_11_gmx_trjconv_str___config
      input_top_path: wf_step_5_grompp/output_tpr_path
      input_structure_path: wf_step_6_mdrun/output_gro_path
    run: biobb/biobb_adapters/cwl/biobb_analysis/gmx_trjconv_str.cwl
    out:
    - output_str_path
cwlVersion: v1.0
class: Workflow
requirements:
  SubworkflowFeatureRequirement:
    class: SubworkflowFeatureRequirement
inputs:
  wf_step_1_setup.yml___setup_step_1_pdb___output_pdb_path: string
  wf_step_1_setup.yml___setup_step_1_pdb___config: string
  wf_step_1_setup.yml___setup_step_4_editconf___config: string
  wf_step_1_setup.yml___setup_step_6_grompp___config: string
  wf_step_1_setup.yml___setup_step_7_genion___config: string
  wf_step_2_min.yml___min_step_1_grompp___config: string
  wf_step_2_min.yml___min_step_1_grompp___output_tpr_path: string
  wf_step_2_min.yml___min_step_3_gmx_energy___config: string
  wf_step_2_min.yml___min_step_3_gmx_energy___output_xvg_path: string
  wf_step_3_nvt.yml___nvt_step_1_grompp___config: string
  wf_step_3_nvt.yml___nvt_step_3_gmx_energy___config: string
  wf_step_3_nvt.yml___nvt_step_3_gmx_energy___output_xvg_path: string
  wf_step_4_npt.yml___npt_step_1_grompp___config: string
  wf_step_4_npt.yml___npt_step_3_gmx_energy___config: string
  wf_step_4_npt.yml___npt_step_3_gmx_energy___output_xvg_path: string
  wf_step_5_grompp___config: string
  wf_step_5_grompp___output_tpr_path: string
  wf_step_6_mdrun___output_gro_path: string
  wf_step_7_gmx_rms___config: string
  wf_step_7_gmx_rms___output_xvg_path: string
  wf_step_8_gmx_rms___config: string
  wf_step_8_gmx_rms___output_xvg_path: string
  wf_step_9_gmx_rgyr___config: string
  wf_step_10_gmx_image___config: string
  wf_step_11_gmx_trjconv_str___config: string
outputs:
  wf_step_1_setup.yml___setup_step_1_pdb___output_pdb_path:
    type: File
    outputSource: wf_step_1_setup.yml/setup_step_1_pdb___output_pdb_path
  wf_step_1_setup.yml___setup_step_2_fix_side_chain___output_pdb_path:
    type: File
    outputSource: wf_step_1_setup.yml/setup_step_2_fix_side_chain___output_pdb_path
  wf_step_1_setup.yml___setup_step_3_pdb2gmx___output_gro_path:
    type: File
    outputSource: wf_step_1_setup.yml/setup_step_3_pdb2gmx___output_gro_path
  wf_step_1_setup.yml___setup_step_3_pdb2gmx___output_top_zip_path:
    type: File
    outputSource: wf_step_1_setup.yml/setup_step_3_pdb2gmx___output_top_zip_path
  wf_step_1_setup.yml___setup_step_4_editconf___output_gro_path:
    type: File
    outputSource: wf_step_1_setup.yml/setup_step_4_editconf___output_gro_path
  wf_step_1_setup.yml___setup_step_5_solvate___output_gro_path:
    type: File
    outputSource: wf_step_1_setup.yml/setup_step_5_solvate___output_gro_path
  wf_step_1_setup.yml___setup_step_5_solvate___output_top_zip_path:
    type: File
    outputSource: wf_step_1_setup.yml/setup_step_5_solvate___output_top_zip_path
  wf_step_1_setup.yml___setup_step_6_grompp___output_tpr_path:
    type: File
    outputSource: wf_step_1_setup.yml/setup_step_6_grompp___output_tpr_path
  wf_step_1_setup.yml___setup_step_7_genion___output_gro_path:
    type: File
    outputSource: wf_step_1_setup.yml/setup_step_7_genion___output_gro_path
  wf_step_1_setup.yml___setup_step_7_genion___output_top_zip_path:
    type: File
    outputSource: wf_step_1_setup.yml/setup_step_7_genion___output_top_zip_path
  wf_step_2_min.yml___min_step_1_grompp___output_tpr_path:
    type: File
    outputSource: wf_step_2_min.yml/min_step_1_grompp___output_tpr_path
  wf_step_2_min.yml___min_step_2_mdrun___output_trr_path:
    type: File
    outputSource: wf_step_2_min.yml/min_step_2_mdrun___output_trr_path
  wf_step_2_min.yml___min_step_2_mdrun___output_gro_path:
    type: File
    outputSource: wf_step_2_min.yml/min_step_2_mdrun___output_gro_path
  wf_step_2_min.yml___min_step_2_mdrun___output_edr_path:
    type: File
    outputSource: wf_step_2_min.yml/min_step_2_mdrun___output_edr_path
  wf_step_2_min.yml___min_step_2_mdrun___output_log_path:
    type: File
    outputSource: wf_step_2_min.yml/min_step_2_mdrun___output_log_path
  wf_step_2_min.yml___min_step_2_mdrun___output_cpt_path:
    type: File
    outputSource: wf_step_2_min.yml/min_step_2_mdrun___output_cpt_path
  wf_step_2_min.yml___min_step_3_gmx_energy___output_xvg_path:
    type: File
    outputSource: wf_step_2_min.yml/min_step_3_gmx_energy___output_xvg_path
  wf_step_3_nvt.yml___nvt_step_1_grompp___output_tpr_path:
    type: File
    outputSource: wf_step_3_nvt.yml/nvt_step_1_grompp___output_tpr_path
  wf_step_3_nvt.yml___nvt_step_2_mdrun___output_trr_path:
    type: File
    outputSource: wf_step_3_nvt.yml/nvt_step_2_mdrun___output_trr_path
  wf_step_3_nvt.yml___nvt_step_2_mdrun___output_gro_path:
    type: File
    outputSource: wf_step_3_nvt.yml/nvt_step_2_mdrun___output_gro_path
  wf_step_3_nvt.yml___nvt_step_2_mdrun___output_edr_path:
    type: File
    outputSource: wf_step_3_nvt.yml/nvt_step_2_mdrun___output_edr_path
  wf_step_3_nvt.yml___nvt_step_2_mdrun___output_log_path:
    type: File
    outputSource: wf_step_3_nvt.yml/nvt_step_2_mdrun___output_log_path
  wf_step_3_nvt.yml___nvt_step_2_mdrun___output_cpt_path:
    type: File
    outputSource: wf_step_3_nvt.yml/nvt_step_2_mdrun___output_cpt_path
  wf_step_3_nvt.yml___nvt_step_3_gmx_energy___output_xvg_path:
    type: File
    outputSource: wf_step_3_nvt.yml/nvt_step_3_gmx_energy___output_xvg_path
  wf_step_4_npt.yml___npt_step_1_grompp___output_tpr_path:
    type: File
    outputSource: wf_step_4_npt.yml/npt_step_1_grompp___output_tpr_path
  wf_step_4_npt.yml___npt_step_2_mdrun___output_trr_path:
    type: File
    outputSource: wf_step_4_npt.yml/npt_step_2_mdrun___output_trr_path
  wf_step_4_npt.yml___npt_step_2_mdrun___output_gro_path:
    type: File
    outputSource: wf_step_4_npt.yml/npt_step_2_mdrun___output_gro_path
  wf_step_4_npt.yml___npt_step_2_mdrun___output_edr_path:
    type: File
    outputSource: wf_step_4_npt.yml/npt_step_2_mdrun___output_edr_path
  wf_step_4_npt.yml___npt_step_2_mdrun___output_log_path:
    type: File
    outputSource: wf_step_4_npt.yml/npt_step_2_mdrun___output_log_path
  wf_step_4_npt.yml___npt_step_2_mdrun___output_cpt_path:
    type: File
    outputSource: wf_step_4_npt.yml/npt_step_2_mdrun___output_cpt_path
  wf_step_4_npt.yml___npt_step_3_gmx_energy___output_xvg_path:
    type: File
    outputSource: wf_step_4_npt.yml/npt_step_3_gmx_energy___output_xvg_path
  wf_step_5_grompp___output_tpr_path:
    type: File
    outputSource: wf_step_5_grompp/output_tpr_path
  wf_step_6_mdrun___output_trr_path:
    type: File
    outputSource: wf_step_6_mdrun/output_trr_path
  wf_step_6_mdrun___output_gro_path:
    type: File
    outputSource: wf_step_6_mdrun/output_gro_path
  wf_step_6_mdrun___output_edr_path:
    type: File
    outputSource: wf_step_6_mdrun/output_edr_path
  wf_step_6_mdrun___output_log_path:
    type: File
    outputSource: wf_step_6_mdrun/output_log_path
  wf_step_6_mdrun___output_cpt_path:
    type: File
    outputSource: wf_step_6_mdrun/output_cpt_path
  wf_step_7_gmx_rms___output_xvg_path:
    type: File
    outputSource: wf_step_7_gmx_rms/output_xvg_path
  wf_step_8_gmx_rms___output_xvg_path:
    type: File
    outputSource: wf_step_8_gmx_rms/output_xvg_path
  wf_step_9_gmx_rgyr___output_xvg_path:
    type: File
    outputSource: wf_step_9_gmx_rgyr/output_xvg_path
  wf_step_10_gmx_image___output_traj_path:
    type: File
    outputSource: wf_step_10_gmx_image/output_traj_path
  wf_step_11_gmx_trjconv_str___output_str_path:
    type: File
    outputSource: wf_step_11_gmx_trjconv_str/output_str_path
