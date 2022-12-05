#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: Download the PDBbind refined database

doc: |-
  Download the PDBbind refined database

baseCommand: python3

hints:
  DockerRequirement:
    dockerPull: ndonyapour/pdbbind_refined_v2020

requirements:
  InlineJavascriptRequirement: {}

inputs:
  script:
    type: string
    inputBinding:
      position: 1
    default: /generate_pdbbind_complex.py

  index_file_name:
    label: The index file name
    type: string
    format:
    - edam:format_2330
    inputBinding:
      prefix: --index_file_name
      position: 2
    default: INDEX_refined_data.2020

  base_dir:
    label: The base_dir path
    type: string
    format:
    - edam:format_2330
    inputBinding:
      prefix: --base_dir
      position: 3
    default: /refined-set

  query:
    label: query str to search the dataset, Pandas query doesn't support slash(/) in column names please use Kd_Ki instead of Kd/Ki
    doc: |-
      query str to search the dataset. Pandas query doesn't support slash(/) in column names please use Kd_Ki instead of Kd/Ki
      Type: string
      File type: input
      Accepted formats: txt
    type: string
    format:
    - edam:format_2330
    inputBinding:
      prefix: --query
      position: 4

  output_txt_path:
    label: Path to the text dataset file
    doc: |-
      Path to the text dataset file
      Type: string
      File type: output
      Accepted formats: txt
    type: string
    format:
    - edam:format_2330
    inputBinding:
      prefix: --output_txt_path
      position: 5
    default: system.log

  output_pdb_paths:
    label: Path to the input file
    doc: |-
      Path to the input file
      Type: string
      File type: input
      Accepted formats: pdb
    type: string
    format:
    - edam:format_1476 # pdb

  output_sdf_paths:
    label: Path to the input file
    doc: |-
      Path to the input file
      Type: string
      File type: input
      Accepted formats: sdf
    type: string
    format:
    - edam:format_3814 # sdf

  min_row:
    label: The row min index
    doc: |-
      The row min inex
      Type: int
    type: int?
    format:
    - edam:format_2330
    inputBinding:
      position: 6
      prefix: --min_row

  max_row:
    label: The row max index
    doc: |-
      The row max inex
      Type: int
    type: int?
    format:
    - edam:format_2330
    inputBinding:
      position: 7
      prefix: --max_row

  convert_Kd_dG:
    label: If this is set to true, dG will be calculated
    doc: If this is set to true, dG will be calculated
    type: string
    format:
    - edam:format_2330
    inputBinding:
      prefix: --convert_Kd_dG
      position: 8
    default: False

  experimental_dGs:
    label: Experimental Free Energies of Binding
    doc: |-
      Experimental Free Energies of Binding
    type: string?
    format:
    - edam:format_2330

outputs:

  output_txt_path:
    label: Path to the txt file
    doc: |-
      Path to the txt file
    type: File
    outputBinding:
      glob: $(inputs.output_txt_path)
    format: edam:format_2330

  output_pdb_paths:
    label: Path to the input file
    doc: |-
      Path to the input file
      Type: string
      File type: input
      Accepted formats: pdb
    type: File[]
    outputBinding:
      glob: ./*.pdb #or  "*.pdb"
    format: edam:format_1476

  output_sdf_paths:
    label: Path to the input file
    doc: |-
      Path to the input file
      Type: string
      File type: input
      Accepted formats: sdf
    type: File[]
    outputBinding:
      glob: ./*.sdf #or  "*.sdf"
    format: edam:format_3814

  experimental_dGs:
    label: Experimental Free Energies of Binding
    doc: |-
      Experimental Free Energies of Binding
    type: ["null", {"type": "array", "items": "float"}]
    outputBinding:
      glob: $(inputs.output_txt_path)
      loadContents: true
      outputEval: |
        ${
          var lines = self[0].contents.split("\n");
          var experimental_dGs = [];
          for (var i = 0; i < lines.length; i++) {
            var words = lines[i].split(" ");
            if (words.length > 2) {
              var experimental_dG = parseFloat(words[2]);
              experimental_dGs.push(experimental_dG);
            }
          }

          if (experimental_dGs.length == 0) {
            return null;
          } else {
            return experimental_dGs;
          }
        }

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl