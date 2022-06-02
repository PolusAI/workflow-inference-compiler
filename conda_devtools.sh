conda install -c conda-forge cwltool graphviz biosimspace

conda install -c conda-forge pytest pytest-cov mypy pylint types-requests types-PyYAML types-setuptools
# NOTE: https://github.com/wearepal/data-science-types has been archived and is
# no longer under active development. So most of the API is covered, but there
# are some functions which are missing stubs.
conda install -c conda-forge data-science-types