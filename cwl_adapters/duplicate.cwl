
#!/usr/bin/env cwl-runner
cwlVersion: v1.1

class: CommandLineTool

label: Duplicates a pdbqt file n times, where n is the length of another array.

doc: |-
  Duplicates a pdbqt file n times, where n is the length of another array.

baseCommand: python3

requirements:
  InlineJavascriptRequirement: {}

inputs:
  input_pdbqt_singleton_path:
    label: Path to the input PDBQT file to be duplicated
    doc: |-
      Path to the input PDBQT file to be duplicated
      Type: string
      File type: input
      Accepted formats: pdbqt
      Example file: https://github.com/bioexcel/biobb_vs/raw/master/biobb_vs/test/data/vina/vina_ligand.pdbqt
    type: File
    format:
    - edam:format_1476

  input_pdbqt_array_path:
    label: Path to the input PDBQT file array
    doc: |-
      Path to the input PDBQT file array
      Type: string
      File type: input
      Accepted formats: pdbqt
      Example file: https://github.com/bioexcel/biobb_vs/raw/master/biobb_vs/test/data/vina/vina_ligand.pdbqt
    type: File[]
    format:
    - edam:format_1476

  output_pdbqt_path:
    label: Path to the output PDBQT files
    doc: |-
      Path to the output PDBQT files
      Type: string
      File type: output
      Accepted formats: pdbqt
      Example file: https://github.com/bioexcel/biobb_vs/raw/master/biobb_vs/test/reference/vina/ref_output_vina.pdbqt
    type: string
    format:
    - edam:format_1476
    default: .

outputs:
  output_pdbqt_path:
    label: Path to the output PDBQT files
    doc: |-
      Path to the output PDBQT files
    type: File[]
    outputBinding:
      outputEval: |
        ${
          var files = [];
          for (var i = 0; i < inputs.input_pdbqt_array_path.length; i++) {
            files.push(inputs.input_pdbqt_singleton_path)
          }
          return files;
        }
    format: edam:format_1476

#stdout: stdout

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
