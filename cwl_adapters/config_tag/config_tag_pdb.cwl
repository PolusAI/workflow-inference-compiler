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
  pdb_id:
    type: string
    format:
    - edam:format_2330

  filter:
    type: string?

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
          config["pdb_code"] = inputs.pdb_id;
          if (inputs.filter) {
            config["filter"] = inputs.filter;
          }
          return JSON.stringify(config);
        }

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
