#!/bin/bash -e
CONDA="Miniforge-pypy3-$(uname)-$(uname -m).sh"
curl -L -O  https://github.com/conda-forge/miniforge/releases/latest/download/"$CONDA"
chmod +x "$CONDA"
./"$CONDA" -b
~/Miniforge-pypy3/bin/mamba init
rm -f "$CONDA"