#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0

label: Returns a dictionary of the given arguments as a JSON-encoded string.
doc: |-
  Returns a dictionary of the given arguments as a JSON-encoded string.

baseCommand: echo # Anything, unused

requirements:
  InlineJavascriptRequirement: {}

inputs:
  water_type:
    type: string
    format:
    - edam:format_2330

  force_field:
    type: string
    format:
    - edam:format_2330

  ignh:
    type: boolean
    format:
    - edam:format_2330

  merge:
    type: boolean
    format:
    - edam:format_2330

# TODO: his

  output_config_string:
    label: A dictionary of the given arguments as a JSON-encoded string.
    doc: |-
      A dictionary of the given arguments as a JSON-encoded string.
    type: string?
    format:
    - edam:format_2330

outputs:
  output_config_string:
    label: A dictionary of the given arguments as a JSON-encoded string.
    doc: |-
      A dictionary of the given arguments as a JSON-encoded string.
    type: string
    #format: edam:format_2330 # "'str' object does not support item assignment""
    outputBinding:
      outputEval: |
        ${
          var config = {};
          config["water_type"] = inputs.water_type;
          config["force_field"] = inputs.force_field;
          config["ignh"] = inputs.ignh;
          config["merge"] = inputs.merge;
          return JSON.stringify(config);
        }

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
