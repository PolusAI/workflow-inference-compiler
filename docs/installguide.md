# Install Guide

## download

Rather than download archive files, it is highly recommended to use the [git](https://git-scm.com) version control system. The software and plugins are currently available via the following http or ssh URLs:

```
https://github.com/PolusAI/workflow-inference-compiler.git
```
```
git@github.com:PolusAI/workflow-inference-compiler.git
```

Regular users can use http, but developers should use ssh.

Either way, the commands

```shell
git clone --recursive <URL>
cd workflow-inference-compiler
```

should automatically download both the software and the plugins. Note that if you just use

```shell
git clone <URL>
cd workflow-inference-compiler
```
it will not download the plugins; you will then need to run the following command:

```
git submodule init && git submodule update
```

Developers who forked the upstream repository should also run the command

```
git remote add upstream https://github.com/PolusAI/workflow-inference-compiler.git
```

## docker

Most of the plugins have been packaged up into [docker](https://www.docker.com) containers. If docker is not installed, you can compile workflows and generate DAGs but the workflows will fail at runtime.

### known issues

When running Docker for Mac, in some cases execution will hang. See [all containers hang](https://github.com/docker/for-mac/issues/5081) and [error waiting for container](https://github.com/docker/for-mac/issues/5139) for details. The current workaround is to simply restart Docker.

If you are experiencing hanging, and if the command `ps aux | grep com.docker | wc -l` returns more than 1000, this is likely the issue. If restarting Docker via the GUI doesn't work, try `sudo pkill com.docker && sudo pkill Docker`.

## conda

[conda](https://en.wikipedia.org/wiki/Conda_(package_manager)) is an open source, cross platform package management system. It can install both Python dependencies and system binary dependencies, so conda is essentially a replacement for `pip`. The associated package distributions named `anaconda` and `miniconda` provide a database of packages that can be used with the `conda` command. (There is also an open source package distribution called [conda-forge](https://conda-forge.org)) Either one works, so if you already have anaconda installed then great. Otherwise, [miniconda](https://docs.conda.io/en/latest/miniconda.html) is all you need.

You can install conda and the system dependencies with the following commands:

```
./install_conda.sh  # install_conda.bat on Windows
conda create --name wic
conda activate wic
./install_system_deps.sh  # install_system_deps.bat on Windows
```

Note that if you close your terminal, the next time you open your terminal you will need to re-activate the wic environment:

```
conda activate wic
```

## install

To install into the wic environment, simply use the following command:

```
pip install -e ".[all]"
```

Developers should also install the git pre-commit hooks:
```
pre-commit install
```

You should now have the `wic` executable available in your terminal.

## Jupyter notebook visualization

The visualization currently needs to be in its own conda environment.

```
./install_conda.sh
conda create --name vis
conda activate vis
./install_nglview.sh
pip install -e ".[all]"
```

## testing

To test your installation, you can run the example in README.md:

```
wic --yaml examples/gromacs/tutorial.yml --run_local --quiet
```

You can also run the automated test suite. This will run the full set of tests, which takes about an hour on a laptop with a GPU.

```
pytest -m serial && pytest -m "not serial" --workers 4
```

If you're in a hurry, just run

```
pytest -m serial
```

which tests the compiler only (not the runtime) and only takes about a minute.

### known issues

At runtime, the following warning message may be printed to the console repeatedly:

```
WARNING: No ICDs were found. Either,
- Install a conda package providing a OpenCL implementation (pocl, oclgrind, intel-compute-runtime, beignet) or
- Make your system-wide implementation visible by installing ocl-icd-system conda package.
```

It does not appear to cause any problems, and moreover the suggested solutions don't seem to work, so it appears that this warning can be ignored.

## documentation

The latest documentation is available on [readthedocs](https://workflow-inference-compiler.readthedocs.io/en/latest/) as well as under the `docs/` folder. To build the documentation in html format, use the commands

```
cd docs/
make html
```

Then simply open the file `docs/_build/html/index.html` using any web browser.

## IntelliSense code completion

You can use any text editor to modify the YAML files, but we *highly* recommend [Visual Studio Code](https://code.visualstudio.com). vscode is an extremely popular development environment written by Microsoft, but it is open source. It has a large number of extensions for pretty much everything. Most importantly, vscode provides [IntelliSense](https://code.visualstudio.com/docs/editor/intellisense) code completion. We now support this for our custom YAML file format. Just simultaneously press ctrl-space (the control key and the spacebar) and a list of possible values will be shown specific to each yaml tag (along with in-line documentation via tooltips). We also support completion for gromacs mdp files!

![code_completion.png](code_completion.png)

To enable this feature, simply install the YAML vscode extension (by Red Hat). Then add the following line to your `.vscode/settings.json` file (or simply use our settings.json file), and run the following command.

```json
"yaml.schemas": {"autogenerated/schemas/wic.json": "*.yml"},
```

```
wic --generate_schemas_only
```

(You should also run this command anytime you change the input/output API. If you are getting validation errors, try re-running this command!)

After a ~10 second delay, vscode should display "Validating against the Workflow Interence Compiler schema" just above the first line in a \*.yml file.

We also *highly* recommend you add the following line to settings.json, which "Controls whether completions should be computed based on words in the document."

```json
"[yaml]": {
    "editor.suggest.showWords": false
},
```

Using the default value of true, IntelliSense will additionally suggest any text whatsoever that you have previously entered in the current yml file (or other yml files). Of course, even if you select a bad suggestion, IntelliSense will immediately validate it against the schema and prominently highlight any invalid text with a red underline. Moreover, the very first thing the compiler does is use a separate validation implementation to check for errors in the yml files. But in this case it's better to simply disable that particular vscode feature.

### known issues

* Strings which are part of enumerated types (e.g. certain gromacs mdp options) are case-sensitive. Unfortunately, there does not seem to be a simple way to allow case-insensitive strings in the schema. This should not be a problem at runtime (we can ignore case in the compiler), but vscode will complain that 'Value must be ...'

### vscode extensions

When you open the project folder, VSCode should prompt you to install the following list of recommended extensions.

These vscode extensions are strongly recommended for users:

* YAML (Red Hat)

These vscode extensions are recommended for developers:

* Python
* MyPy
* Remote Development
* Docker
* CWL (Rabix/Benten)
* Git Extension Pack
* Github Actions (Mathieu Dutour)
* autoDocstring
* Protein Viewer