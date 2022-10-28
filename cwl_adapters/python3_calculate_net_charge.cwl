#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: Calculate the Total Charge 

doc: |-
  Calculate the total charge of a given ligand

baseCommand: python3
#arguments: ["../scripts/examples/calculate_net_charge.py"]

# hints:
#   DockerRequirement:
#     dockerPull: jakefennick/scripts
requirements:
  InlineJavascriptRequirement: {}

inputs:
  script:
    type: File
    format:
    - edam:format_2330
    inputBinding:
      position: 1
      #default: ../scripts/examples/calculate_net_charge.py
  input_path:
    label: Path to the input file
    doc: |-
      Path to the input file
      Type: string
      File type: input
      Accepted formats: pdb, mol2
    type: File
    format:
    - edam:format_1476 # pdb
    - edam:format_3816 # mol2
    inputBinding:
      position: 2
      prefix: --input_path

  addhydrogens:
    label: adds hydrogens to the system 
    doc: adds hydrogens to the system 
    type: boolean?
    format:
    - edam:format_2330
    inputBinding:
      prefix: --addhydrogens
      position: 3
    default: False

  output_log_path:
    label: Path to the log file
    doc: |-
      Path to the log file
      Type: string
      File type: output
      Accepted formats: log
    type: string
    format:
    - edam:format_2330
    inputBinding:
      prefix: --output_log_path
      position: 4
    default: system.log

  net_charge:
    label: Calculated total charge
    doc: |-
      Calculated total charge
    type: string
    format:
    - edam:format_2330

outputs:
  output_log_path:
    label: Path to the log file
    doc: |-
      Path to the log file
    type: File
    outputBinding:
      glob: $(inputs.output_log_path)
    format: edam:format_2330

  net_charge:
    label: Calculated total charge
    doc: |-
      Calculated total charge
    type: int
    outputBinding:
      glob: $(inputs.output_log_path)
      loadContents: true
      outputEval: |
        ${
          var lines = self[0].contents.split("\n");
          // The correct line should be of the form
          // Calculated total charge: 1
          var net_charge_line = lines.filter(function(s) {return s.split(" ")[0] == "Calculated"})[0];
          var net_charge_string = net_charge_line.split(" ").filter(function(s) {return !isNaN(parseFloat(s))})[0];
          var net_charge = parseInt(net_charge_string);
          return net_charge
        }


$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl