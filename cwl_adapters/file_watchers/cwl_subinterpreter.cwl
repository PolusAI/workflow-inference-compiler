#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: Speculatively executes an arbitrary CommandLineTool (upto max_times) by file watching / polling --cachedir. This is primarily intended for parsing logfiles before the associated CWL process has finished.

doc: |-
  Speculatively executes an arbitrary CommandLineTool (upto max_times) by file watching / polling --cachedir. This is primarily intended for parsing logfiles before the associated CWL process has finished.

baseCommand: cwl_subinterpreter

inputs:
  cachedir_path:
    label: The full absolute path to the --cachedir cwltool directory. You should also use this same directory when invoking RealtimePlots.py
    doc: |
      The full absolute path to the --cachedir cwltool directory. You should also use this same directory when invoking RealtimePlots.py
    type: string
    format:
    - edam:format_2330 # 'Textual format'
    inputBinding:
      position: 1
      prefix: --cachedir_path

  file_pattern:
    label: Filenames that match this pattern will be watched / polled for changes. i.e. '*.log'
    doc: |
      Filenames that match this pattern will be watched / polled for changes. i.e. '*.log'
    type: string
    format:
    - edam:format_2330 # 'Textual format'
    inputBinding:
      position: 2
      prefix: --file_pattern

  cwl_tool:
    label: The filename (without .cwl extension) of the CommandLineTool to speculatively execute.
    doc: |
      The filename (without .cwl extension) of the CommandLineTool to speculatively execute.
    type: string
    format:
    - edam:format_2330 # 'Textual format'
    inputBinding:
      position: 3
      prefix: --cwl_tool

  max_times:
    label: The maximum number of times to speculatively execute cwl_tool. This is used to guarantee termination in case of errors.
    doc: |
      The maximum number of times to speculatively execute cwl_tool. This is used to guarantee termination in case of errors.
    type: string
    format:
    - edam:format_2330 # 'Textual format'
    inputBinding:
      position: 4
      prefix: --max_times

  config:
    label: A JSON-encoded string which contains the arguments (i.e. the `in:` tag) to the wrapped CommandLineTool. Make sure to escape any substrings as necessary!
    doc: |-
      A JSON-encoded string which contains the arguments (i.e. the `in:` tag) to the wrapped CommandLineTool. Make sure to escape any substrings as necessary!
    type: string
    format:
    - edam:format_2330 # 'Textual format'
    inputBinding:
      position: 5
      prefix: --config

  root_workflow_yml_path:
    label: The full absolute path to the root workflow yml file.
    doc: |
      The full absolute path to the root workflow yml file.
    type: string
    format:
    - edam:format_2330 # 'Textual format'
    inputBinding:
      position: 6
      prefix: --root_workflow_yml_path

  homedir:
    label: The full absolute path to the users home directory.
    doc: |
      The full absolute path to the root users home directory.
    type: string
    format:
    - edam:format_2330 # 'Textual format'
    inputBinding:
      position: 7
      prefix: --homedir

outputs:
  output_log_path:
    label: Path to the output log file
    doc: |-
      Path to the output log file
    type: File
    outputBinding:
      glob: $(inputs.cwl_tool)_only.log
    format: edam:format_2330 # 'Textual format'

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
