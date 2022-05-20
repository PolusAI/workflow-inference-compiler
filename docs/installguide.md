## download

Rather than download archive files, it is highly recommended to use the [git](https://git-scm.com) version control system. The software and plugins are currently available via the following http or ssh URLs:

```
https://github.com/jfennick/workflow_inference_compiler.git

git@github.com:jfennick/workflow_inference_compiler.git
```

Regular users can use http, but developers should use ssh.

Either way, the command `git clone --recursive <URL>` should automatically download both the software and the plugins.

Note that if you just use `git clone <URL>` it will not download the plugins; you will then need to run the following command:

```
git submodule init && git submodule update
```

## docker

Most of the plugins have been packaged up into [docker](https://www.docker.com) containers. If docker is not installed, you can compile workflows and generate DAGs but the workflows will fail at runtime.

## conda

[conda](https://en.wikipedia.org/wiki/Conda_(package_manager)) is an open source, cross platform package management system. It can install both Python dependencies and system binary dependencies, so conda is essentially a replacement for `pip`. The associated package distributions named `anaconda` and `miniconda` provide a database of packages that can be used with the `conda` command. (There is also an open source package distribution called [conda-forge](https://conda-forge.org)) Either one works, so if you already have anaconda installed then great. Otherwise, [miniconda](https://docs.conda.io/en/latest/miniconda.html) is all you need.

After conda is installed, enter the following commands:

```
conda create --name wic
conda activate wic
./conda_devtools.sh
```

Note that if you close your terminal, the next time you open your terminal you will need to re-activate the wic environment:

```
conda activate wic
```

## install

To install into the wic environment, simply use the following command:

```
pip install .
```

You should now have the `wic` executable available in your terminal.

## testing

To test your installation, you can run the example in the [readme](../README.md#quick-start):

```
wic --yaml examples/gromacs/tutorial.yml --cwl_dir . --cwl_run_local True
```

You can also run the automated test suite. This will run the full set of tests, which takes about 20-30 minutes on a laptop.

```
pytest
```

If you're in a hurry, just run

```
pytest -m fast
```

which tests the compiler only (not the runtime) and only takes about 15 seconds.

## vscode

You can use any text editor to modify the YAML files, but we highly recommend [Visual Studio Code](https://code.visualstudio.com). vscode is an extremely popular development environment written by Microsoft, but it is open source. It has a large number of extensions for pretty much everything. Most importantly, vscode provides [IntelliSense](https://code.visualstudio.com/docs/editor/intellisense) code completion. In the very near future, we plan on creating an IntelliSense extension to support our custom YAML format.

### extensions

I recommend the following vscode extensions for developers:

* Python
* MyPy
* Remote Development
* Docker
* CWL (Rabix/Benten)