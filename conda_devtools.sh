conda install -y -c conda-forge -c schrodinger pymol-bundle
# If you want to use the GUI, also install
# pip install PyQt5
# At runtime, CWL uses the Docker image jakefennick/scripts

# Comment out biosimspace because it is a massive dependency,
# and for the very limited use case of file format conversions, we don't use
# enough of the API to justify installing it for IDE support.
#conda install -y -c conda-forge -c michellab biosimspace
# At runtime, CWL uses the Docker image jakefennick/biosimspace

conda install -y -c conda-forge cwltool nodejs graphviz openbabel mdanalysis
# NOTE: cwltool needs nodejs for InlineJavascriptRequirement

conda install -y -c conda-forge pytest pytest-cov pytest-parallel mypy pylint types-requests types-PyYAML types-setuptools
# NOTE: https://github.com/wearepal/data-science-types has been archived and is
# no longer under active development. So most of the API is covered, but there
# are some functions which are missing stubs.
conda install -y -c conda-forge data-science-types

conda install -y -c conda-forge wget

# Temporarily use conda (instead of the Docker image jakefennick/autodock_vina)
# due to issues with some AWS EC2 instances.
conda install -y -c conda-forge vina
