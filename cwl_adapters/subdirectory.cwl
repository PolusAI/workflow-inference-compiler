#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: Accesses a subdirectory within a directory

doc: |
  Accesses a subdirectory within a directory

baseCommand: "true"

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
    - $(inputs.directory)

inputs:
  directory:
    type: Directory

  glob_pattern:
    type: string
    default: "."

outputs:
  subdirectory:
    type: Directory
    outputBinding:
      glob: $(inputs.glob_pattern)

stdout: stdout

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
