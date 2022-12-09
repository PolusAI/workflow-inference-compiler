cwlVersion: v1.0

class: CommandLineTool

label: Run a Bash script

doc: |
  Run a Bash script

baseCommand: bash
requirements:
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
    - $(inputs.script)
    - $(inputs.docker_file)
    - $(inputs.python_script)

inputs:
  script:
    type: File
    format:
    - edam:format_2330
    inputBinding:
      position: 1

  docker_file:
    type: File
    format:
    - edam:format_2330
    inputBinding:
      position: 2

  tag:
    type: string
    format:
    - edam:format_3816
    inputBinding:
      position: 3

  python_script:
    type: File
    format:
    - edam:format_2330
    inputBinding:
      position: 4

outputs:
  stdout:
    type: File
    format: edam:format_2330
    outputBinding:
      glob: stdout

stdout: stdout

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
