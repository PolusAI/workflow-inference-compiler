# Comment out pymol-bundle because it conflicts with `pip install toil[cwl]` below.
#conda install -y -c conda-forge -c schrodinger pymol-bundle
# If you want to use the GUI, also install
# pip install PyQt5
# At runtime, CWL uses the Docker image jakefennick/scripts

# Comment out biosimspace because it is a massive dependency,
# and for the very limited use case of file format conversions, we don't use
# enough of the API to justify installing it for IDE support.
#conda install -y -c conda-forge -c michellab biosimspace
# At runtime, CWL uses the Docker image jakefennick/biosimspace

conda install -y -c conda-forge nodejs graphviz openbabel mdanalysis
# NOTE: cwltool needs nodejs for InlineJavascriptRequirement

# "Warning: Could not load "/miniconda/bin/../lib/graphviz/libgvplugin_pango.so.6"
#  - It was found, so perhaps one of its dependents was not.  Try ldd."
# See https://github.com/conda-forge/graphviz-feedstock/issues/35#issuecomment-786368065
conda install -y -c conda-forge xorg-libxrender

# NOTE: https://github.com/wearepal/data-science-types has been archived and is
# no longer under active development. So most of the API is covered, but there
# are some functions which are missing stubs.
conda install -y -c conda-forge data-science-types

conda install -y -c conda-forge wget
conda install -y -c conda-forge zip # Not installed by default on ubuntu

# NOTE: The [cwl] extra installs an embedded cwltool within toil-cwl-runner.
# You can NOT `conda install cwltool` and then `pip install toil` !
conda install -y -c conda-forge pip
pip install 'toil[cwl]'

# The ruptures library needs to compile its binary wheel during pip install
# Even though the compilers are already installed
# (i.e. the below command prints "All requested packages are already installed.")
# if you don't explicitly install the compilers package,
# the ruptures binary wheel build will fail with
# "#  include <stdlib.h>"
# ...
# "ERROR: Failed building wheel for ruptures"
conda install -y -c conda-forge compilers
