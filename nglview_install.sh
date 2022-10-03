# nglview has a transitive dependency on mistune==2.0.4 which conflicts with
# cwltool's transitive dependency of mistune==0.8.4.
# Fortunately, due to the file watching strategy in src/vis/viewer.py does not
# need to know the structure of the workflows, nglview and cwltool do not
# need to be installed into the same conda environment.

conda install -y -c conda-forge nglview

# mdtraj needs binary build dependencies (specifically cython) so do not install via pip
conda install -y -c conda-forge mdtraj

conda install -y pip

pip install nglview jupyterlab ipywidgets==7.7.1
# See https://github.com/nglviewer/nglview/issues/1032

jupyter-nbextension enable nglview --py --sys-prefix