#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: Globs filenames from disk into a CWL File metadata array. Useful for scattering in a later step. Does not load any file contents.

doc: |
  Globs filenames from disk into a CWL File metadata array. Useful for scattering in a later step. Does not load any file contents.

baseCommand: echo # Anything, unused

requirements:
  InlineJavascriptRequirement: {}

inputs:
  directory:
    type: string
    inputBinding:
      position: 1

  glob_pattern:
    type: string
    inputBinding:
      position: 2

  output_file_format:
    type: string

outputs:
  output_files_path:
    type: File[]
    outputBinding:
      glob: $(inputs.directory)$(inputs.glob_pattern)
      outputEval: |
        ${
          var files = [];
          for (var i = 0; i < self.length; i++) {
            self[i]["format"] = inputs.output_file_format;
            files.push(self[i]);
          }
          return files;
        }

  output_filenames_path:
    type: string[]
    outputBinding:
      glob: $(inputs.directory)$(inputs.glob_pattern)
      outputEval: |
        ${
          var filenames = [];
          for (var i = 0; i < self.length; i++) {
            filenames.push(self[i].basename); //basename returns the filename (i.e. with extension)
          }
          return filenames;
        }

stdout: stdout

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
