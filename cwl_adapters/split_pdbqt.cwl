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
arguments: [$(inputs.input_path), "-o", "pdbqt", "-O", $(inputs.output_pdb_path), "-m", "-xh", "-xn"]
# NOTE: These arguments must be given individually; they cannot be concatenated together.
# (e.g. -xrhn) Otherwise, all but the first argument will be silently ignored!

# -m Produces multiple output files, to allow:
#     Splitting: e.g.        obabel infile.mol -O new.smi -m
#       puts each molecule into new1.smi new2.smi etc
# ...
# pdbqt  AutoDock PDBQT format
# Reads and writes AutoDock PDBQT (Protein Data Bank, Partial Charge (Q), & Atom Type (T)) format
# Note that the torsion tree is by default. Use the ``r`` write option
# to prevent this.

#Read Options, e.g. -ab
#  b  Disable automatic bonding
#  d  Input file is in dlg (AutoDock docking log) format

#Write Options, e.g. -xr
#  b  Enable automatic bonding
#  r  Output as a rigid molecule (i.e. no branches or torsion tree)
#  c  Combine separate molecular pieces of input into a single rigid molecule (requires "r" option or will have no effect)
#  s  Output as a flexible residue
#  p  Preserve atom indices from input file (default is to renumber atoms sequentially)
#  h  Preserve hydrogens
#  n  Preserve atom names

# NOTE: The version of openbabel in this container is old; This may or may not cause problems.
# See https://github.com/openbabel/openbabel/issues/2435
hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/biobb_chemistry:3.7.0--pyhdfd78af_0

requirements:
  InlineJavascriptRequirement: {}

inputs:
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

  output_pdb_path:
    label: Path to the output file
    doc: |-
      Path to the output file
      Type: string?
      File type: output
      Accepted formats: pdb
    type: string?
    format:
    - edam:format_1476 # pdb
    default: system.pdbqt

outputs:
  output_pdb_path:
    label: Path to the output file
    doc: |-
      Path to the output file
    type: File[]
    outputBinding:
      glob: "$(inputs.output_pdb_path.slice(0, -6))*.pdbqt" # e.g. "ligand.pdb" -> "ligand*.pdb"
    format: edam:format_1476 # pdb

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
