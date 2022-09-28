#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: This class is a wrapper of the Open Babel tool.

doc: |-
  Small molecule format conversion for structures or trajectories. Open Babel is a chemical toolbox designed to speak the many languages of chemical data. It's an open, collaborative project allowing anyone to search, convert, analyze, or store data from molecular modeling, chemistry, solid-state materials, biochemistry, or related areas. Visit the official page.

baseCommand: obabel
#Usage: 
#obabel[-i<input-type>] <infilename> [-o<output-type>] -O<outfilename> [Options]
#...
#Options, other than -i -o -O -m, must come after the input files.
arguments: [$(inputs.input_path), "-o", "sdf", "-O", $(inputs.output_sdf_path), "-r"]
# NOTE: These arguments must be given individually; they cannot be concatenated together.
# (e.g. -xrhn) Otherwise, all but the first argument will be silently ignored!

# "-r	Remove all but the largest contiguous fragment (strip salts)"

#sdf  MDL MOL format
#Reads and writes V2000 and V3000 versions

#Open Babel supports an extension to the MOL file standard
#that allows cis/trans and tetrahedral stereochemistry to be
#stored in 0D MOL files. The tetrahedral stereochemistry is
#stored as the atom parity, while the cis/trans stereochemistry
#is stored using Up and Down bonds similar to how it is
#represented in a SMILES string. Use the ``S`` option
#when reading or writing if you want to avoid storing
#or interpreting stereochemistry in 0D MOL files.

#Read Options, e.g. -as
# s  determine chirality from atom parity flags
#       The default setting for 2D and 3D is to ignore atom parity and
#       work out the chirality based on the bond
#       stereochemistry (2D) or coordinates (3D).
#       For 0D the default is already to determine the chirality
#       from the atom parity.
# S  do not read stereochemistry from 0D MOL files
#       Open Babel supports reading and writing cis/trans
#       and tetrahedral stereochemistry to 0D MOL files.
#       This is an extension to the standard which you can
#       turn off using this option.
# T  read title only
# P  read title and properties only
#       When filtering an sdf file on title or properties
#       only, avoid lengthy chemical interpretation by
#       using the ``T`` or ``P`` option together with the
#       :ref:`copy format <Copy_raw_text>`.

#Write Options, e.g. -x3
# 3  output V3000 not V2000 (used for >999 atoms/bonds) 
# a  write atomclass if available
# m  write no properties
# w  use wedge and hash bonds from input (2D only)
# v  always specify the valence in the valence field
#      The default behavior is to only specify the valence if it
#      is not consistent with the MDL valence model.
#      So, for CH4 we don't specify it, but we do for CH3.
#      This option may be useful to preserve the correct number of
#      implicit hydrogens if a downstream tool does not correctly
#      implement the MDL valence model (but does honor the valence
#      field).
# S  do not store cis/trans stereochemistry in 0D MOL files
# A  output in Alias form, e.g. Ph, if present
# E  add an ASCII depiction of the molecule as a property
# H  use HYD extension (always on if mol contains zero-order bonds)

hints:
  DockerRequirement:
    #dockerPull: jakefennick/scripts
    dockerPull: quay.io/biocontainers/biobb_chemistry:3.7.0--pyhdfd78af_0

inputs:
  first_molecule:
    label: Index of the first molecule (1-based)
    doc: |-
      Input Index of the first molecule (1-based)
      Type: string
    type: string?
    format:
    - edam:format_2330 # 'Textual format'
    inputBinding:
      prefix: -f
    default: 1

  last_molecule:
    label: Index of the last molecule (1-based)
    doc: |-
      Input Index of the last molecule (1-based)
      Type: string
    type: string?
    format:
    - edam:format_2330 # 'Textual format'
    inputBinding:
      prefix: -l
    default: 1

  input_path:
    label: Path to the input file
    doc: |-
      Path to the input file
      Type: string
      File type: input
      Accepted formats: dat, ent, fa, fasta, gro, inp, log, mcif, mdl, mmcif, mol, mol2, pdb, pdbqt, png, sdf, smi, smiles, txt, xml, xtc
      Example file: https://github.com/bioexcel/biobb_chemistry/raw/master/biobb_chemistry/test/data/babel/babel.smi
    type: File
    format:
    - edam:format_1637
    - edam:format_1476
    - edam:format_1929
    - edam:format_1929
    - edam:format_2033
    - edam:format_3878
    - edam:format_2030
    - edam:format_1477
    - edam:format_3815
    - edam:format_1477
    - edam:format_3815
    - edam:format_3816
    - edam:format_1476
    - edam:format_1476
    - edam:format_3603
    - edam:format_3814
    - edam:format_1196
    - edam:format_1196
    - edam:format_2033
    - edam:format_2332
    - edam:format_3875

  output_sdf_path:
    label: Path to the output file
    doc: |-
      Path to the output file
      Type: string
      File type: output
      Accepted formats: sdf
    type: string
    format:
    - edam:format_3814 # sdf
    default: system.sdf

outputs:
  output_sdf_path:
    label: Path to the output file
    doc: |-
      Path to the output file
    type: File
    outputBinding:
      glob: $(inputs.output_sdf_path)
    format: edam:format_3814 # sdf

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
