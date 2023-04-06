#!/bin/bash -e
# NOTE: mamba is a drop-in replacement for conda, just much faster.
# (i.e. You can replace mamba with conda below.)
# See https://github.com/conda-forge/miniforge#mambaforge-pypy3
CONDA=conda
if [ "$(which mamba)" ]; then
    CONDA=mamba
fi

# mdtraj needs binary build dependencies (specifically cython) so do not install via pip
$CONDA install -y -c conda-forge mdtraj

# mdtraj and pypy CANNOT currently be installed into the same environment:
# $CONDA install -y -c conda-forge mdtraj pypy

$CONDA install -y pip

pip install jupyterlab ipytree

# See release notes about jupyterlab
# https://github.com/nglviewer/nglview/releases/tag/v3.0.4
pip install nglview

jupyter-nbextension enable nglview --py --sys-prefix
