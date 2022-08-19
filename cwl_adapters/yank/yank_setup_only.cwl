#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: Calculate binding free energy of receptor ligand systems using Yank.

doc: |-
  Calculate binding free energy of receptor ligand systems using Yank.

requirements:
  - class: DockerRequirement
    dockerPull: docker.io/jakefennick/yank
# See https://www.commonwl.org/user_guide/15-staging/index.html
# The documentation on this could be better. It is indeed true that
# by default, every input file is staged into its own read-only temporary
# directory. Enabling InitialWorkDirRequirement will stage all of the input
# files into the output directory, which is writeable. Note that even if you
# do NOT enable InitialWorkDirRequirement, you can still write out files into
# subdirectories of the temporary working directory /private/tmp/docker_tmp*
# However, it appears that only files that are in the root working directory are
# copied to the output directory; subdirectories (and any files contained) are not!

#  --move-outputs        Move output files to the workflow output directory and delete intermediate output directories (default).
#  --leave-outputs       Leave output files in intermediate output directories.
#  --copy-outputs        Copy output files to the workflow output directory and don't delete intermediate output directories.

# (Notice how the documentation says files, but doesn't mention subdirectories.)
# "Did not find output file with glob pattern: '['output/setup/setup.log']'."
# Thus, these files will NOT be accessible to glob!!! (glob looks for files
# w.r.t. the output directory.) Perhaps $(runtime.tmpdir) can be used, but thus
# far I have been unable to change the default value of {}. The solution is to
# enable InitialWorkDirRequirement, which eliminates any need to copy, i.e.
# it is a workaround for the failure to copy files in subdirectories.
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing: |
      ${
        var lst = [];
        for (var i = 0; i < inputs.input_dir_path.length; i++) {
          var dict = {
            "entry": inputs.input_dir_path[i],
            "writable": true // Important!
          };
          lst.push(dict);
        }
        //lst.push(inputs.yaml);
        return lst;
      }
# NOTE: Yank uses the following snippet to determine the location of the checkpoint files:
# with moltools.utils.temporary_cd(self._script_dir):
#            self._check_resume()
# In other words, in addition to staging the checkpoint files, we also need to
# stage the yaml script file so that Yank looks for the checkpoint files in the
# SAME directory!!!

baseCommand: python3
arguments: [$(inputs.script)]
inputs:
  script:
    type: File
    default:
      class: File
      location: yank_wrapper.py

  input_dir_path:
    label: Output directory
    type:
      type: array
      items: [Directory, File]
    format: edam:format_2330 # 'Textual format'
    # NOTE: Cannot provide default value here; add [] to the yml file.
    #default: []
    # No inputBinding

  output_dir_path:
    label: Output directory
    type: string
    format: edam:format_2330 # 'Textual format'
    default: output
    # No inputBinding

  phase:
    type: string
    format: edam:format_string
    default: production
    inputBinding:
      prefix: --phase

  yaml:
    label: Input YAML script
    type: File
    format: edam:format_3750
    inputBinding:
      prefix: --yaml

# NOTE: It is unclear if cwltool supports mutually exclusive parameters.
# See https://github.com/common-workflow-language/cwltool/issues/358
# See https://www.commonwl.org/user_guide/misc/#setting-mutually-exclusive-parameters
# For now, we can create a separate CommandLineTool for the receptor/ligand setup route.

#  input_receptor_path:
#    label: Input receptor pdb file
#    type: File
#    format:
#    - edam:format_1476
#    inputBinding:
#      prefix: --input_receptor_path

#  input_ligand_path:
#    label: Input ligand mol2 file
#    type: File
#    format:
#    - edam:format_3816
#    inputBinding:
#      prefix: --input_ligand_path

# See http://getyank.org/latest/yamlpages/systems.html#yaml-systems-user-defined

  input_complex_top_zip_path:
    label: Input complex top zip file
    type: File
    format:
    - edam:format_3987
    inputBinding:
      prefix: --input_complex_top_zip_path

  input_complex_crd_path:
    label: Input ligand gro file
    type: File
    format:
    - edam:format_2033
    inputBinding:
      prefix: --input_complex_crd_path

  input_ligand_top_zip_path:
    label: Input ligand top zip file
    type: File
    format:
    - edam:format_3987
    inputBinding:
      prefix: --input_ligand_top_zip_path

  input_ligand_crd_path:
    label: Input ligand gro file
    type: File
    format:
    - edam:format_2033
    inputBinding:
      prefix: --input_ligand_crd_path

outputs:
  output_log_path:
    label: Output log file
    type: File
    format: edam:format_2330 # 'Textual format'
    outputBinding:
      glob: output/setup/setup.log # $(runtime.outdir)/output/setup/setup.log also works

# Actually, these outputs are only relevant for the receptor/ligand setup route.

#  output_complex_setup_crd_path:
#    label: Output coordinates file for the complex phase (AMBER crd)
#    type: File
#    format: edam:format_3878
#    outputBinding:
#      glob: output/setup/systems/**/complex.inpcrd

#  output_complex_setup_top_path:
#    label: Output topology file for the complex phase (AMBER crd)
#    type: File
#    format: edam:format_3881
#    outputBinding:
#      glob: output/setup/systems/**/complex.prmtop

#  output_solvent_setup_crd_path:
#    label: Output coordinates file for the solvent phase (AMBER crd)
#    type: File
#    format: edam:format_3878
#    outputBinding:
#      glob: output/setup/systems/**/solvent.inpcrd

#  output_solvent_setup_top_path:
#    label: Output topology file for the solvent phase (AMBER crd)
#    type: File
#    format: edam:format_3881
#    outputBinding:
#      glob: output/setup/systems/**/solvent.prmtop

# See https://rabix.io/cwl-patterns.html
# Note that we need to list individual output files above (so we can add format
# & refer to them in later workflow steps), but by default CWL will copy them
# from their subdirectory to the root output directory. To maintain
# subdirectories, we also need to use output_dir_path!
  output_dir_path:
    type:
      type: array
      items: [Directory, File]
    outputBinding: {glob: "*"}
    # When InitialWorkDirRequirement is enabled, the entire docker_tmp*
    # working directory gets copied to the output directory.
    # Use "*" to glob into it, or use "." to output everything.
    format: edam:format_2330 # 'Textual format'

#  output_dir_path:
#    type: Directory
#    format: edam:format_2330 # 'Textual format'
#    outputBinding: { glob: "*" }
     # This version of output_dir_path works with ".", but does not work with "*":
     # "Multiple matches for output item that is a single file."

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl