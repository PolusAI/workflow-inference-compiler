:: NOTE: mamba is a drop-in replacement for conda, just much faster.
:: (i.e. You can replace mamba with conda below.)
:: See https://github.com/conda-forge/miniforge#mambaforge-pypy3
set CONDA=conda
where /q mamba
if not ERRORLEVEL 1 (set CONDA=mamba)

:: %CONDA% clean --all --yes
%CONDA% env update --file system_deps_windows.yml
:: pip cache purge