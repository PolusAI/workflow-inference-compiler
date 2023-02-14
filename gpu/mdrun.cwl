#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: Wrapper of the GROMACS mdrun module with GPUs.

doc: |-
  MDRun is the main computational chemistry engine within GROMACS. It performs Molecular Dynamics simulations, but it can also perform Stochastic Dynamics, Energy Minimization, test particle insertion or (re)calculation of energies.

baseCommand: gmx
arguments: ["mdrun"]

hints:
  cwltool:CUDARequirement:
    cudaVersionMin: "11.4"
    cudaComputeCapabilityMin: "3.0"
    cudaDeviceCount: 0
  DockerRequirement:
    dockerPull: gromacs/gromacs:2021.4 # uses CUDA 10.1 which is compatible with our EC2 instances.
    # NOTE: None of the nvcr.io/hpc/gromacs images seem to work

inputs:

  nb_terms:
    label: Device to perform non-bonded interactions on
    doc: |-
      Device to perform non-bonded interaction
    type: string?
    inputBinding:
      prefix: -nb
      position: 1
    default: gpu

  bonded_terms:
    label: Device to perform bonded interactions on
    doc: |-
      Device to perform bonded interaction
    type: string?
    inputBinding:
      prefix: -bonded
      position: 2
    default: gpu

  pme_terms:
    label: Device to perform PME calculations on
    doc: |-
      Device to perform PME calculations on
    type: string?
    inputBinding:
      prefix: -pme
      position: 3
    default: gpu

  number_threads_openmp:
    label: Number of openmp threads to start
    doc: |-
      Number of openmp threads to start
    type: int
# Cannot set total number of threads -nt with the gromacs/gromacs images, else:
# "Fatal error:
# Setting the total number of threads is only supported with thread-MPI and
# GROMACS was compiled without thread-MPI"
    inputBinding:
      prefix: -ntomp
      position: 4
    default: 2 # We want to essentially disable the CPU-GPU load-balancing.
    # However, since some operations are CPU only, using only 1 CPU will
    # slow down overall runtime (see Amdahl's Law). However, we also do NOT
    # want to just use all CPUs, because if we run many simulations in
    # parallel this will massively oversubscribe the CPUs.
    # N^2 oversubscription is not okay; 2*N oversubscription is okay.

  input_tpr_path:
    label: Path to the portable binary run input file TPR
    doc: |-
      Path to the portable binary run input file TPR
      Type: string
      File type: input
      Accepted formats: tpr
      Example file: https://github.com/bioexcel/biobb_md/raw/master/biobb_md/test/data/gromacs/mdrun.tpr
    type: File
    format:
    - edam:format_2333
    inputBinding:
      position: 5
      prefix: -s

  output_trr_path:
    label: Path to the GROMACS uncompressed raw trajectory file TRR
    doc: |-
      Path to the GROMACS uncompressed raw trajectory file TRR
      Type: string
      File type: output
      Accepted formats: trr
      Example file: https://github.com/bioexcel/biobb_md/raw/master/biobb_md/test/reference/gromacs/ref_mdrun.trr
    type: string
    format:
    - edam:format_3910
    inputBinding:
      position: 6
      prefix: -o
    default: system.trr

  output_crd_path:
    label: Path to the output GROMACS structure GRO file
    doc: |-
      Path to the output GROMACS structure GRO file
      Type: string
      File type: output
      Accepted formats: gro
      Example file: https://github.com/bioexcel/biobb_md/raw/master/biobb_md/test/reference/gromacs/ref_mdrun.gro
    type: string
    format:
    - edam:format_2033
    inputBinding:
      position: 7
      prefix: -c
    default: system.gro

  output_edr_path:
    label: Path to the output GROMACS portable energy file EDR
    doc: |-
      Path to the output GROMACS portable energy file EDR
      Type: string
      File type: output
      Accepted formats: edr
      Example file: https://github.com/bioexcel/biobb_md/raw/master/biobb_md/test/reference/gromacs/ref_mdrun.edr
    type: string
    format:
    - edam:format_2330
    inputBinding:
      position: 8
      prefix: -e
    default: system.edr

  output_log_path:
    label: Path to the output GROMACS trajectory log file LOG
    doc: |-
      Path to the output GROMACS trajectory log file LOG
      Type: string
      File type: output
      Accepted formats: log
      Example file: https://github.com/bioexcel/biobb_md/raw/master/biobb_md/test/reference/gromacs/ref_mdrun.log
    type: string
    format:
    - edam:format_2330
    inputBinding:
      position: 9
      prefix: -g
    default: system.log

  input_cpt_path:
    label: Path to the input GROMACS checkpoint file CPT
    doc: |-
      Path to the input GROMACS checkpoint file CPT
      Type: string
      File type: input
      Accepted formats: cpt
      Example file: null
    type: File?
    format:
    - edam:format_2333
    inputBinding:
      prefix: -cpi

  output_xtc_path:
    label: Path to the GROMACS compressed trajectory file XTC
    doc: |-
      Path to the GROMACS compressed trajectory file XTC
      Type: string
      File type: output
      Accepted formats: xtc
      Example file: null
    type: string
    format:
    - edam:format_3875
    inputBinding:
      prefix: -x
    default: system.xtc

  output_cpt_path:
    label: Path to the output GROMACS checkpoint file CPT
    doc: |-
      Path to the output GROMACS checkpoint file CPT
      Type: string
      File type: output
      Accepted formats: cpt
      Example file: null
    type: string
    format:
    - edam:format_2333
    inputBinding:
      prefix: -cpo
    default: system.cpt

  output_dhdl_path:
    label: Path to the output dhdl.xvg file only used when free energy calculation
      is turned on
    doc: |-
      Path to the output dhdl.xvg file only used when free energy calculation is turned on
      Type: string
      File type: output
      Accepted formats: xvg
      Example file: null
    type: string
    format:
    - edam:format_2033
    inputBinding:
      prefix: -dhdl
    default: system.xvg

outputs:
  output_trr_path:
    label: Path to the GROMACS uncompressed raw trajectory file TRR
    doc: |-
      Path to the GROMACS uncompressed raw trajectory file TRR
    type: File
    outputBinding:
      glob: $(inputs.output_trr_path)
    format: edam:format_3910

  output_crd_path:
    label: Path to the output GROMACS structure GRO file
    doc: |-
      Path to the output GROMACS structure GRO file
    type: File
    outputBinding:
      glob: $(inputs.output_crd_path)
    format: edam:format_2033

  output_edr_path:
    label: Path to the output GROMACS portable energy file EDR
    doc: |-
      Path to the output GROMACS portable energy file EDR
    type: File
    outputBinding:
      glob: $(inputs.output_edr_path)
    format: edam:format_2330

  output_log_path:
    label: Path to the output GROMACS trajectory log file LOG
    doc: |-
      Path to the output GROMACS trajectory log file LOG
    type: File
    outputBinding:
      glob: $(inputs.output_log_path)
    format: edam:format_2330

  output_xtc_path:
    label: Path to the GROMACS compressed trajectory file XTC
    doc: |-
      Path to the GROMACS compressed trajectory file XTC
    type: File?
    outputBinding:
      glob: $(inputs.output_xtc_path)
    format: edam:format_3875

  output_cpt_path:
    label: Path to the output GROMACS checkpoint file CPT
    doc: |-
      Path to the output GROMACS checkpoint file CPT
    type: File?
    outputBinding:
      glob: $(inputs.output_cpt_path)
    format: edam:format_2333

  output_dhdl_path:
    label: Path to the output dhdl.xvg file only used when free energy calculation
      is turned on
    doc: |-
      Path to the output dhdl.xvg file only used when free energy calculation is turned on
    type: File?
    outputBinding:
      glob: $(inputs.output_dhdl_path)
    format: edam:format_2033

$namespaces:
  edam: https://edamontology.org/
  cwltool: http://commonwl.org/cwltool#

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
