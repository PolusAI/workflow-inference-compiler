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

# mdtraj and pypy CAN be installed into the same environment, but they need to be
# installed first and simultaneously, so the dependency solver can find a solution:
# $CONDA install -y -c conda-forge mdtraj pypy
# The solver currently finds pypy3.6 and mdtraj 1.9.4, i.e. installing into the
# same environment as wic will downgrade us from python 3.9 to 3.6, and we
# cannot simply call both install scripts sequentially.
# But we do not want to require installing molecular modeling dependencies
# into wic, so we should just use two separate environments anyway.

$CONDA install -y pip

pip install jupyterlab ipytree

# See https://github.com/nglviewer/nglview/pull/1045
# Since nglview appears to be abandoned, I cherry-picked these two commits so
# we don't have to depend on an external github user for ipywidgets 8 support.
pip install "git+https://github.com/jfennick/nglview.git" ipywidgets

# Alternatively, we can pin to the previous version, but this will generate
# dependency conflicts over time.
# pip install nglview ipywidgets==7.7.1

jupyter-nbextension enable nglview --py --sys-prefix
