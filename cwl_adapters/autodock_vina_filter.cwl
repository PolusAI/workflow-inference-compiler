
#!/usr/bin/env cwl-runner
cwlVersion: v1.1

class: CommandLineTool

label: Filters results of the AutoDock Vina software.

doc: |-
  This class applies a cutoff to the docking scores of the AutoDock Vina software.

baseCommand: python3

hints:
  DockerRequirement:
    dockerPull: jakefennick/scripts

requirements:
  InlineJavascriptRequirement: {}

inputs:
  script:
    type: string
    inputBinding:
      position: 1
    default: /autodock_vina_filter.py # NOTE: Initial / required

# NOTE: To make inference work, at least one of the following two log inputs needs to be non-optional.
# However, since there is no way in CWL (I think) to make inputs mutually exclusive, we would need to
# supply fake default values in the form of anonymous Files and then perform a check in the script.
# However, there is (at least one) bug in cwltool related to anonymous files (files that start with _:...) + Docker, etc.

  input_log_path:
    label: Path to the log file
    doc: |-
      Path to the log file
      Type: string
      File type: output
      Accepted formats: log
      Example file: https://github.com/bioexcel/biobb_vs/raw/master/biobb_vs/test/reference/vina/ref_output_vina.log
    type: File?
    format:
    - edam:format_2330
    inputBinding:
      position: 2
      prefix: --input_log_path
    #default: {class: File, basename: nonexistent_logfile.log, contents: "", format: edam:format_2330}

  input_log_paths:
    label: Path to the log files
    doc: |-
      Path to the log files
      Type: string
      File type: output
      Accepted formats: log
      Example file: https://github.com/bioexcel/biobb_vs/raw/master/biobb_vs/test/reference/vina/ref_output_vina.log
    type: File[]?
    format:
    - edam:format_2330
    inputBinding:
      position: 2
      prefix: --input_log_paths
    #default: [{class: File, basename: nonexistent_logfile.log, contents: "", format: edam:format_2330}]

  docking_score_cutoff:
    label: Cutoff threshold for filtering docking scores
    doc: |-
      Cutoff threshold for filtering docking scores
      Type: float
    type: float
    format:
    - edam:format_2330
    inputBinding:
      position: 3
      prefix: --docking_score_cutoff

  max_num_decoys_per_ligand:
    label: Maximum number of decoys per initial ligand
    doc: |-
      Maximum number of decoys per initial ligand
      Type: int
    type: int
    format:
    - edam:format_2330
    inputBinding:
      position: 4
      prefix: --max_num_decoys_per_ligand
    default: -1

  max_num_decoys_total:
    label: Maximum number of decoys total
    doc: |-
      Maximum number of decoys total
      Type: int
    type: int
    format:
    - edam:format_2330
    inputBinding:
      position: 5
      prefix: --max_num_decoys_total
    default: -1

  input_batch_pdbqt_path:
    label: Path to the input PDBQT ligands
    doc: |-
      Path to the input PDBQT ligands
      Type: string
      File type: input
      Accepted formats: pdbqt
      Example file: https://github.com/bioexcel/biobb_vs/raw/master/biobb_vs/test/data/vina/vina_ligand.pdbqt
    type: {"type": "array", "items": {"type": "array", "items": "File"}}
    #type: File[][] # Invalid syntax; [] syntactic sugar only works for 1D arrays.
    format:
    - edam:format_1476
#    inputBinding:
#      position: 6 # Since the type is File[], this means starting at this position.
#      prefix: --input_batch_pdbqt_path

# NOTE: This is only used so we can create explicit edges.
# The scatter-related inference bugs are now sorted out, so this can probably be removed.
  output_batch_pdbqt_path:
    label: Path to the output PDBQT batch files
    doc: |-
      Path to the output PDBQT batch files
      Type: string
      File type: output
      Accepted formats: pdbqt
      Example file: https://github.com/bioexcel/biobb_vs/raw/master/biobb_vs/test/reference/vina/ref_output_vina.pdbqt
    type: string
    format:
    - edam:format_1476
#    inputBinding:
#      position: 7
#      prefix: --output_batch_pdbqt_path
    default: .

outputs:
# NOTE: If docking_score_cutoff is too negative and filters out all of the files,
# you will get the following runtime error message:
# ("Error collecting output for parameter 'output_batch_pdbqt_path': cwl_adapters/autodock_vina_filter.cwl:73:3: 'NoneType' object does not support item assignment", {})
  output_batch_pdbqt_path:
    label: Path to the output PDBQT file
    doc: |-
      Path to the output PDBQT file
# While 2D array inputs seem to be working, all of my attempts to output a 2D array using
# outputEval have failed due to generating a stack trace in the cwltool python process:
# ("Error collecting output for parameter 'output_batch_pdbqt_path': cwl_adapters/autodock_vina_filter.cwl:134:3: list indices must be integers or slices, not str", {})
# I have even tried outputEval: $([[]]) and outputEval: $(inputs.output_batch_pdbqt_path)
# Note that outputEval: $([]) and outputEval: $(inputs.output_batch_pdbqt_path[0]) with type: File[] works.
    #type: {"type": "array", "items": {"type": "array", "items": "File"}}
    #type: File[][] # See above
    type: File[]
    outputBinding:
      glob: indices.txt # This determines what binds to self[0]
      loadContents: true # If true, this additionally binds self[0].contents
      # NOTE: According to the CWL specification, loadContents only reads the first 64KB.
      # In version 1.2, any more is an error; in versions 1.0 and 1.1, the file contents is silently truncated!
      # Since log files and/or pdb files can easily exceed that, we do all
      # processing using an external script, read in the final indices here,
      # and index into the existing json. This preserves the additional json metadata fields,
      # (location, basename, class, checksum, size, format, path)
      # which would be difficult to recover from an external script which is only given path.
      # (An external script can write a cwl.output.json file, which will alternatively determine the outputs.)
      outputEval: |
        ${
          var lines = self[0].contents.split("\n");
          var files = []; // In this case, flatten the 2D nested array into a 1D array
          for (var i = 0; i < lines.length; i++) {
            var indices = lines[i].split(" ");
            //var docking_score = parseFloat(indices[0]);
            var mol_idx = parseInt(indices[1]);
            var mode_idx = parseInt(indices[2]);
            var file = inputs.input_batch_pdbqt_path[mol_idx][mode_idx];
            files.push(file);
          }
          return files;
        }
    format: edam:format_1476

stdout: stdout

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
