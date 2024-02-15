# Workflow Inference Compiler

![buid-status](https://readthedocs.org/projects/workflow-inference-compiler/badge/?version=latest&style=svg)

Scientific computing can be difficult in practice due to various complex software issues. In particular, chaining together software packages into a computational pipeline can be very error prone. Using the [Common Workflow Language](https://www.commonwl.org) (CWL) greatly helps, but like many other workflow languages users still need to explicitly specify how to connect inputs & outputs. The Workflow Inference Compiler allows users to specify computational protocols at a very high level of abstraction, it automatically infers almost all connections between inputs & outputs, and it compiles to CWL for execution.

![BPS Poster](BPS_poster.svg)

## Documentation
The documentation is available on [readthedocs](https://workflow-inference-compiler.readthedocs.io/en/latest/).
## Example Workflows
The main examples are chosen from classical molecular dynamics but like CWL, the compiler is general purpose and is not limited to any specific domain.

[Molecular Modeling Workflows](https://github.com/PolusAI/mm-workflows)
## Quick Start
See the [installation guide](docs/installguide.md) for more details, but:
```
git clone https://github.com/PolusAI/workflow-inference-compiler.git
cd workflow-inference-compiler

cd install/
./install_conda.sh  # install_conda.bat on Windows
source ~/.bashrc  # skip on Windows
conda create --name wic
conda activate wic
./install_pypy.sh  # Optional, skip on Arm Macs
./install_system_deps.sh  # install_system_deps.bat on Windows
cd ..

pip install -e ".[all]"
pre-commit install  # Required for developers

cd install
./install_biobb_adapters.sh
./install_mm-workflows.sh
./install_image-workflows.sh
cd ..
```
```
wic --generate_schemas_only
wic --yaml ../mm-workflows/examples/gromacs/tutorial.yml --run_local --quiet
```
That last command will infer edges, compile the yml to CWL, generate a GraphViz diagram of the root workflow, and run it locally.

```yaml
label: Conjugate Gradient
steps:
  - grompp:
      in:
        config:
          mdp:
            integrator: cg
            nsteps: 1000
  - mdrun:
      in:
        # Use GPU by default
        bonded_terms: cpu
        pme_terms: cpu
  - gmx_energy:
      in:
        config:
          terms: [Potential]
        output_xvg_path: energy_min_cg.xvg
```
The subworkflow [`examples/gromacs/cg.yml`](https://github.com/PolusAI/mm-workflows/blob/main/examples/gromacs/cg.yml) in `mm-workflows` is shown above, and the GraphViz diagram of the root workflow [`examples/gromacs/tutorial.yml`](https://github.com/PolusAI/mm-workflows/blob/main/examples/gromacs/tutorial.yml) in `mm-workflows` is shown below.

![Workflow](examples/gromacs/tutorial.yml.gv.png)

If you add the --parallel flag to the above command then, in another terminal, you can view the plots in real-time:
```
conda activate wic
cd install && ./install_timeseriesplots.sh && cd ..
timeseriesplots
```

![Plots](examples/gromacs/plots.png)

You can also view the 3D structures in the Jupyter notebook `src/vis/viewer.ipynb`. The visualization currently needs to be in its own conda environment.

```
install/install_conda.sh
conda create --name vis
conda activate vis
install/install_nglview.sh
pip install -e ".[all]"
```

![Plots](docs/tree_viewer.png)

The Workflow Inference Compiler is a [Domain Specific Language](https://en.wikipedia.org/wiki/Domain-specific_language) (DSL) based on the [Common Workflow Language](https://www.commonwl.org). CWL is fantastic, but explicitly constructing the Directed Acyclic Graph (DAG) associated with a non-trivial workflow is not so simple. For example, the workflow shown above is based on the following [gromacs tutorial](https://mmb.irbbarcelona.org/biobb/availability/tutorials/cwl). Instead of writing raw CWL, you can write your workflows in a much simpler yml DSL. For technical reasons edge inference is far from unique, so ***`users should always check that edge inference actually produces the intended DAG`***.

## Edge Inference

The key feature is that in most cases, you do not need to specify any of the edges! They will be automatically inferred for you based on types, file formats, and naming conventions. For more information, see the [user guide](docs/userguide.md#edge-inference-algorithm) If for some reason edge inference fails, there is a syntax for creating [explicit edges](docs/userguide.md#explicit-edges).

## Subworkflows

Subworkflows are very useful for creating reusable, composable building blocks. As shown above, recursive subworkflows are fully supported, and the edge inference algorithm has been very carefully constructed to work across subworkflow boundaries.

## Explicit CWL

Since the yml DSL files are automatically compiled to CWL, users should not have to know any CWL. However, the yml DSL is secretly CWL that is simply missing almost all of the tags! In other words, the compiler merely adds missing information to the files, and so if you know CWL you are free to explicitly add the information yourself. Thus, the yml DSL is intentionally a [leaky abstraction](https://en.wikipedia.org/wiki/Leaky_abstraction).
