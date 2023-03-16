#!/bin/bash -e
# NOTE: mamba is a drop-in replacement for conda, just much faster.
# (i.e. You can replace mamba with conda below.)
# See https://github.com/conda-forge/miniforge#mambaforge-pypy3
CONDA=conda
if [ "$(which mamba)" ]; then
    CONDA=mamba
fi

#$CONDA clean --all --yes

# pypy is ~2X faster than the regular python interpreter.
# We need to install it first so the dependency solver installs it bundled with python 3.9
# (pypy is not yet compatible with 3.10 and 3.11)
$CONDA install -y -c conda-forge pypy "python<3.10"

# Comment out pymol-bundle because it conflicts with `pip install toil[cwl]` below.
#mamba install -y -c conda-forge -c schrodinger pymol-bundle
# If you want to use the GUI, also install
# pip install PyQt5
# At runtime, CWL uses the Docker image jakefennick/scripts

# Comment out biosimspace because it is a massive dependency,
# and for the very limited use case of file format conversions, we don't use
# enough of the API to justify installing it for IDE support.
#$CONDA install -y -c conda-forge -c michellab biosimspace
# At runtime, CWL uses the Docker image jakefennick/biosimspace

$CONDA install -y -c conda-forge nodejs graphviz # openbabel mdanalysis
# NOTE: cwltool needs nodejs for InlineJavascriptRequirement

# "Warning: Could not load "/miniconda/bin/../lib/graphviz/libgvplugin_pango.so.6"
#  - It was found, so perhaps one of its dependents was not.  Try ldd."
# See https://github.com/conda-forge/graphviz-feedstock/issues/35#issuecomment-786368065
$CONDA install -y -c conda-forge xorg-libxrender

# NOTE: https://github.com/wearepal/data-science-types has been archived and is
# no longer under active development. So most of the API is covered, but there
# are some functions which are missing stubs.
$CONDA install -y -c conda-forge data-science-types

$CONDA install -y -c conda-forge wget
$CONDA install -y -c conda-forge zip # Not installed by default on ubuntu

$CONDA install -y -c conda-forge shellcheck

# NOTE: The [cwl] extra installs an embedded cwltool within toil-cwl-runner.
# You can NOT `conda install cwltool` and then `pip install toil` !
# NOTE: ruamel.yaml (a dependency of toil/cwtool) requires compilers
$CONDA install -y -c conda-forge compilers
$CONDA install -y -c conda-forge pip
#pip cache purge