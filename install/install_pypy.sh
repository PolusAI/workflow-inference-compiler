#!/bin/bash -e
# NOTE: mamba is a drop-in replacement for conda, just much faster.
# (i.e. You can replace mamba with conda below.)
# See https://github.com/conda-forge/miniforge#Miniforge-pypy3
CONDA=conda
if [ "$(which mamba)" ]; then
    CONDA=mamba
fi

#$CONDA clean --all --yes
$CONDA env update --file pypy.yml
#pip cache purge