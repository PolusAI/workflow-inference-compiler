#!/usr/bin/env cwl-runner
steps:
  prod_step_1_grompp:
    in:
      config: prod_step_1_grompp___config
      output_tpr_path: prod_step_1_grompp___output_tpr_path
      input_gro_path: prod_step_1_grompp___input_gro_path
      input_top_zip_path: prod_step_1_grompp___input_top_zip_path
    run: biobb/biobb_adapters/cwl/biobb_md/grompp.cwl
    out:
    - output_tpr_path
  prod_step_2_mdrun:
    in:
      output_gro_path: prod_step_2_mdrun___output_gro_path
      input_tpr_path: prod_step_1_grompp/output_tpr_path
    run: biobb/biobb_adapters/cwl/biobb_md/mdrun.cwl
    out:
    - output_trr_path
    - output_gro_path
    - output_edr_path
    - output_log_path
    - output_xtc_path
    - output_cpt_path
    - output_dhdl_path
cwlVersion: v1.0
class: Workflow
inputs:
  prod_step_1_grompp___config: string
  prod_step_1_grompp___output_tpr_path: string
  prod_step_1_grompp___input_gro_path: File
  prod_step_1_grompp___input_top_zip_path: File
  prod_step_2_mdrun___output_gro_path: string
outputs:
  prod_step_1_grompp___output_tpr_path:
    type: File
    outputSource: prod_step_1_grompp/output_tpr_path
  prod_step_2_mdrun___output_trr_path:
    type: File
    outputSource: prod_step_2_mdrun/output_trr_path
  prod_step_2_mdrun___output_gro_path:
    type: File
    outputSource: prod_step_2_mdrun/output_gro_path
  prod_step_2_mdrun___output_edr_path:
    type: File
    outputSource: prod_step_2_mdrun/output_edr_path
  prod_step_2_mdrun___output_log_path:
    type: File
    outputSource: prod_step_2_mdrun/output_log_path
  prod_step_2_mdrun___output_cpt_path:
    type: File
    outputSource: prod_step_2_mdrun/output_cpt_path
