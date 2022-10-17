
#!/usr/bin/env cwl-runner
cwlVersion: v1.1

# NOTE: This file is nearly identical to autodock_vina_batch with the primary difference that
# --input_batch_pdbqt_path is replaced with --input_ligand_pdbqt_path.
# (For no obvious reason, --score_only only works with --ligand, not --batch)

class: CommandLineTool

label: Wrapper of the AutoDock Vina software.

doc: |-
  This class performs docking of the ligand to a set of grids describing the target protein via the AutoDock Vina software.

baseCommand: vina # NOTE: Only version >=1.2 supports --batch!
arguments:
- "--autobox"
# NOTE: The documentation for --score_only claims "search space can be omitted" which is not quite correct;
# if you omit --center_* and --size_* you get "ERROR: Grid box dimensions must be greater than 0 Angstrom."
# However, adding --autobox works.

hints:
  DockerRequirement:
    dockerPull: jakefennick/autodock_vina

requirements:
  InlineJavascriptRequirement: {}

inputs:
  input_ligand_pdbqt_path:
    label: Path to the input PDBQT ligand
    doc: |-
      Path to the input PDBQT ligand
      Type: string
      File type: input
      Accepted formats: pdbqt
      Example file: https://github.com/bioexcel/biobb_vs/raw/master/biobb_vs/test/data/vina/vina_ligand.pdbqt
    type: File
    format:
    - edam:format_1476
    inputBinding:
      position: 1
      prefix: --ligand

  input_receptor_pdbqt_path:
    label: Path to the input PDBQT receptor
    doc: |-
      Path to the input PDBQT receptor
      Type: string
      File type: input
      Accepted formats: pdbqt
      Example file: https://github.com/bioexcel/biobb_vs/raw/master/biobb_vs/test/data/vina/vina_receptor.pdbqt
    type: File
    format:
    - edam:format_1476
    inputBinding:
      position: 2
      prefix: --receptor

  local_only:
    label: Do local search only
    doc: Do local search only
    type: boolean?
    #format:
    #- edam_format_2330 # textual format
    inputBinding:
      prefix: --local_only

  score_only:
    label: Do not do any conformational search; simply rescore.
    doc: Do not do any conformational search; simply rescore.
    type: boolean?
    #format:
    #- edam_format_2330 # textual format
    inputBinding:
      prefix: --score_only

# NOTE: This is only used so we can create explicit edges.
# The scatter-related inference bugs are now sorted out, so this can probably be removed.
  #output_batch_dir_path:
  output_batch_pdbqt_path:
    label: Path to the output PDBQT batch directory
    doc: |-
      Path to the output PDBQT batch directory
      Type: string
      File type: output
      Accepted formats: pdbqt
      Example file: https://github.com/bioexcel/biobb_vs/raw/master/biobb_vs/test/reference/vina/ref_output_vina.pdbqt
    type: string?
    format:
    - edam:format_1476
#    inputBinding:
#      position: 4
#      prefix: --dir
#    default: .

  output_log_path:
    label: Path to the log file
    doc: |-
      Path to the log file
      Type: string
      File type: output
      Accepted formats: log
      Example file: https://github.com/bioexcel/biobb_vs/raw/master/biobb_vs/test/reference/vina/ref_output_vina.log
    type: string
    format:
    - edam:format_2330
    default: system.log

  docking_score:
    label: Estimated Free Energy of Binding
    doc: |-
      Estimated Free Energy of Binding
    type: string
    format:
    - edam:format_2330

outputs:
  output_batch_pdbqt_path:
    label: Path to the output PDBQT files
    doc: |-
      Path to the output PDBQT files
    #type: File[]
    type:
      type: array
      items: File
    outputBinding:
      #glob: $(inputs.output_batch_dir_path)/*.pdbqt
      glob: ./*.pdbqt # Use ./* because leading *'s are reserved syntax for Yaml aliases.
    format: edam:format_1476

  output_log_path:
    label: Path to the log file
    doc: |-
      Path to the log file
    type: File
    outputBinding:
      glob: $(inputs.output_log_path)
    format: edam:format_2330

  docking_score:
    label: Estimated Free Energy of Binding
    doc: |-
      Estimated Free Energy of Binding
    type: float
    outputBinding:
      glob: $(inputs.output_log_path)
      loadContents: true
      outputEval: |
        ${
          var lines = self[0].contents.split("\n");
          // The correct line should be of the form
          // Estimated Free Energy of Binding   : -6.053 (kcal/mol) [=(1)+(2)+(3)+(4)]
          var bfe_line = lines.filter(function(s) {return s.split(" ")[0] == "Estimated"})[0];
          var docking_score_string = bfe_line.split(" ").filter(function(s) {return !isNaN(parseFloat(s))})[0];
          var docking_score = parseFloat(docking_score_string);
          return docking_score
        }

stdout: $(inputs.output_log_path)

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
