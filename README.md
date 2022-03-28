# Workflow Inference Compiler

## Quick Start
```
pip install .
wic --yaml examples/gromacs/tutorial.yml --cwl_dir . --cwl_output_intermediate_files True --cwl_run_local True
```

![Workflow](examples/gromacs/tutorial.yml.gv.png)

The Workflow Inference Compiler is a [Domain Specific Language](https://en.wikipedia.org/wiki/Domain-specific_language) (DSL) based on the [Common Workflow Language](https://www.commonwl.org). CWL is fantastic, but explicitly constructing the Directed Acyclic Graph (DAG) associated with a non-trivial workflow is not so simple. For example, the workflow shown above is based on the following [gromacs tutorial](http://mmb.irbbarcelona.org/webdev/slim/biobb/public/availability/tutorials/cwl). Instead of writing raw CWL, you can write your workflows in a much simpler yml DSL.

## Edge Inference

The key feature is that in most cases, you do not need to specify any of the edges! They will be automatically inferred for you based on file formats and other information. The above command will infer edges, compile the yml to CWL, generate a GraphViz diagram of the workflow, and run it locally.

### Maths

Every DAG has a [topological ordering](https://en.wikipedia.org/wiki/Topological_sorting). Since CWL workflows are DAGs, there must be an associated topological ordering. However, since the input yml DSL only contains a linear sequence of steps and does not contain any edge information, we merely have a [linear ordering](https://en.wikipedia.org/wiki/Total_order). The challenge is to promote the linear ordering to a topological ordering by inferring all of the edges. Since we are initially missing information this is far from unique, so ***`users should always check that edge inference actually produces the intended DAG`***.

## Explicit Edges

If for some reason edge inference fails, you can always explicitly specify the edges using `'&var'` and `'*var'` notation. Simply use `'&var'` to create a reference to an output filename and then, in an input in any later step, use `'*var'` to dereference the filename and create an explicit edge between the output and the input. See examples/gromacs for a concrete example. Unfortunately, due to yml syntax, these variables will need to be in quotes.

## Subworkflows

Subworkflows are very useful for creating reusable, composable building blocks. As shown above, recursive subworkflows are fully supported, and the edge inference algorithm has been very carefully constructed to work across subworkflow boundaries. If there are subworkflows, the linear order in the parent workflow is determined by inlineing the subworkflows. (If you are unsure, simply enable `--cwl_inline_subworkflows`) Note that since CWL files for subworkflows are also automatically generated, any subworkflow can be treated as a black box if desired.

### Embedding Independence
Moreover, in addition to ignoring irrelevant details inside a subworkflow, users should also be able to ignore the details of the parent workflow in which a subworkflow is embedded. In other words, subworkflows should be context-free and have 'embedding independence'. Again, the edge inference algorithm has been carefully constructed such that the edges inferred within a subworkflow do not depend on its embedding within a parent workflow. (This requirement is guaranteed by the regression tests!)

## Explicit CWL

Since the yml DSL files are automatically compiled to CWL, users should not have to know any CWL. However, the yml DSL is secretly CWL that is simply missing almost all of the tags! In other words, the compiler merely adds missing information to the files, and so if you know CWL you are free to explicitly add the information yourself. Thus, the yml DSL is intentionally a [leaky abstraction](https://en.wikipedia.org/wiki/Leaky_abstraction).