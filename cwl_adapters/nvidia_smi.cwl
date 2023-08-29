cwlVersion: v1.0

class: CommandLineTool

baseCommand: nvidia-smi

hints:
  cwltool:CUDARequirement:
    cudaVersionMin: "11.4"
    cudaComputeCapability: "3.0"
    cudaDeviceCountMin: 1
    cudaDeviceCountMax: 1
  DockerRequirement:
    dockerPull: nvidia/cuda:11.4.3-base-ubuntu20.04

requirements:
  InlineJavascriptRequirement: {}

inputs:
  query:
    label: display listed info
    doc: |-
      display listed info
      Type: string
      File type: input
      Accepted formats: txt
    type: string
    format:
    - edam:format_2330
    inputBinding:
      separate: false
      prefix: --query-gpu=
      position: 1

  file_format:
    label: output file format
    doc: |-
      output file format
      Type: string
      File type: input
      Accepted formats: txt
    type: string
    format:
    - edam:format_2330
    inputBinding:
      separate: false
      prefix: --format=
      position: 2

  output_txt_path:
    label: the name of ouputfile
    doc: |-
      the name of ouputfile
      Type: string
      File type: input
      Accepted formats: txt
    type: string
    format:
    - edam:format_2330
    inputBinding:
      separate: false
      prefix: --filename=
      position: 3

outputs:
  output_txt_path:
    label: Path to the txt file
    doc: |-
      Path to the txt file
    type: File
    outputBinding:
      glob: $(inputs.output_txt_path)
    format: edam:format_2330

$namespaces:
  edam: https://edamontology.org/
  cwltool: http://commonwl.org/cwltool#

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
