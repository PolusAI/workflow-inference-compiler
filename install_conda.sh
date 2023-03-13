#!/bin/bash -e
CONDA="Mambaforge-pypy3-$(uname)-$(uname -m).sh"
curl -L -O  https://github.com/conda-forge/miniforge/releases/latest/download/"$CONDA"
chmod +x "$CONDA"
./"$CONDA" -b
rm -f "$CONDA"