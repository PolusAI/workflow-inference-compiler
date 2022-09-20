
#!/usr/bin/env cwl-runner
cwlVersion: v1.1

class: CommandLineTool

label: Wrapper of the AutoDock Vina software.

doc: |-
  This class performs docking of the ligand to a set of grids describing the target protein via the AutoDock Vina software.

baseCommand: vina # NOTE: Only version >=1.2 supports --batch!
arguments:
# Need to parse box.pdb and pass in each number separately.
# REMARK BOX CENTER:     0.102     0.019    -0.004 SIZE:    30.195    31.940    27.005
- "--dir" # Need to explicitly pass --dir . in --batch mode
- "."
- "--center_x"
- $(inputs.input_box_path.contents.split("\n")[0].split(" ").filter(function(s) {return !isNaN(parseFloat(s))})[0])
- "--center_y"
- $(inputs.input_box_path.contents.split("\n")[0].split(" ").filter(function(s) {return !isNaN(parseFloat(s))})[1])
- "--center_z"
- $(inputs.input_box_path.contents.split("\n")[0].split(" ").filter(function(s) {return !isNaN(parseFloat(s))})[2])
- "--size_x"
- $(inputs.input_box_path.contents.split("\n")[0].split(" ").filter(function(s) {return !isNaN(parseFloat(s))})[3])
- "--size_y"
- $(inputs.input_box_path.contents.split("\n")[0].split(" ").filter(function(s) {return !isNaN(parseFloat(s))})[4])
- "--size_z"
- $(inputs.input_box_path.contents.split("\n")[0].split(" ").filter(function(s) {return !isNaN(parseFloat(s))})[5])
# NOTE: Cannot use a single javascript expression to create the entire arguments list because CWL treats it as a string:
# "the `arguments` field is not valid because value is a str"
#  ${
#    var words = inputs.input_box_path.contents.split("\n")[0].split(" ");
#    var nums = words.filter(function(s) {return !isNaN(parseFloat(s))});
#    var args = {};
#    args.push("--dir"); // Need to explicitly pass --dir . in --batch mode
#    args.push(".");
#    args.push("--center_x");
#    args.push(nums[0]);
#    args.push("--center_y");
#    args.push(nums[1]);
#    args.push("--center_z");
#    args.push(nums[2]);
#    args.push("--size_x");
#    args.push(nums[3]);
#    args.push("--size_y");
#    args.push(nums[4]);
#    args.push("--size_z");
#    args.push(nums[5]);
#    return args;
#  }

hints:
  DockerRequirement:
    dockerPull: jakefennick/autodock_vina

requirements:
  InlineJavascriptRequirement: {}

inputs:
  input_batch_pdbqt_path:
    label: Path to the input PDBQT ligands
    doc: |-
      Path to the input PDBQT ligands
      Type: string
      File type: input
      Accepted formats: pdbqt
      Example file: https://github.com/bioexcel/biobb_vs/raw/master/biobb_vs/test/data/vina/vina_ligand.pdbqt
    type: File[]
    format:
    - edam:format_1476
    inputBinding:
      position: 1
      prefix: --batch

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

  input_box_path:
    label: Path to the PDB containig the residues belonging to the binding site
    doc: |-
      Path to the PDB containig the residues belonging to the binding site
      Type: string
      File type: input
      Accepted formats: pdb
      Example file: https://github.com/bioexcel/biobb_vs/raw/master/biobb_vs/test/data/vina/vina_box.pdb
    type: File
    format:
    - edam:format_1476
#    inputBinding:
#      position: 3
#      prefix: --input_box_path
    loadContents: true


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

stdout: $(inputs.output_log_path)

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
