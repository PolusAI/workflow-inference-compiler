#!/usr/bin/env cwl-runner
steps:
  npt_step_1_grompp:
    in:
      config: npt_step_1_grompp___config
      input_gro_path: npt_step_1_grompp___input_gro_path
      input_top_zip_path: npt_step_1_grompp___input_top_zip_path
    run: biobb/biobb_adapters/cwl/biobb_md/grompp.cwl
    out:
    - output_tpr_path
  npt_step_2_mdrun:
    run: biobb/biobb_adapters/cwl/biobb_md/mdrun.cwl
    in:
      input_tpr_path: npt_step_1_grompp/output_tpr_path
    out:
    - output_trr_path
    - output_gro_path
    - output_edr_path
    - output_log_path
    - output_xtc_path
    - output_cpt_path
    - output_dhdl_path
  npt_step_3_gmx_energy:
    in:
      config: npt_step_3_gmx_energy___config
      output_xvg_path: npt_step_3_gmx_energy___output_xvg_path
      input_energy_path: npt_step_2_mdrun/output_edr_path
    run: biobb/biobb_adapters/cwl/biobb_analysis/gmx_energy.cwl
    out:
    - output_xvg_path
cwlVersion: v1.0
class: Workflow
inputs:
  npt_step_1_grompp___config: string
  npt_step_1_grompp___input_gro_path: File
  npt_step_1_grompp___input_top_zip_path: File
  npt_step_3_gmx_energy___config: string
  npt_step_3_gmx_energy___output_xvg_path: string
outputs:
  npt_step_1_grompp___output_tpr_path:
    type: File
    outputSource: npt_step_1_grompp/output_tpr_path
  npt_step_2_mdrun___output_trr_path:
    type: File
    outputSource: npt_step_2_mdrun/output_trr_path
  npt_step_2_mdrun___output_gro_path:
    type: File
    outputSource: npt_step_2_mdrun/output_gro_path
  npt_step_2_mdrun___output_edr_path:
    type: File
    outputSource: npt_step_2_mdrun/output_edr_path
  npt_step_2_mdrun___output_log_path:
    type: File
    outputSource: npt_step_2_mdrun/output_log_path
  npt_step_2_mdrun___output_cpt_path:
    type: File
    outputSource: npt_step_2_mdrun/output_cpt_path
  npt_step_3_gmx_energy___output_xvg_path:
    type: File
    outputSource: npt_step_3_gmx_energy/output_xvg_path
