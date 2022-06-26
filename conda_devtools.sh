conda install -c conda-forge wget libarchive # bsdtar is installed with libarchive; see download_data.sh

conda install -c conda-forge -c schrodinger pymol-bundle
# If you want to use the GUI, also install
# pip install PyQt5

conda install -c conda-forge -c michellab biosimspace

conda install -c conda-forge cwltool graphviz openbabel mdanalysis

conda install -c conda-forge pytest pytest-cov mypy pylint types-requests types-PyYAML types-setuptools
# NOTE: https://github.com/wearepal/data-science-types has been archived and is
# no longer under active development. So most of the API is covered, but there
# are some functions which are missing stubs.
conda install -c conda-forge data-science-types