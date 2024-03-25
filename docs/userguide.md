# User Guide

Recall for a moment the vague instructions your PhD advisor hastily scribbled onto the chalkboard about how to do a calculation. Now imagine that those scribbles were actually executable! That's the goal! The goal is to allow high-level, domain-specific concepts to be directly specified in a user-friendly YAML format. Then the abstract scientific protocol is automatically translated into specific concrete steps, executed on a remote job cluster, and automated analysis is performed.

## Main Features / Design Overview

See [overview](overview.md)

## auto-discovery

Many software packages have a way of automatically discovering files which they can use. (examples: [pytest](https://docs.pytest.org/en/latest/explanation/goodpractices.html#conventions-for-python-test-discovery) [pylint](https://pylint.pycqa.org/en/latest/user_guide/usage/run.html))

By default, wic will recursively search for tools / workflows within the directories (and subdirectories) listed in the config file's json tags `search_paths_cwl` and `search_paths_yml`. The paths listed can be absolute or relative. The default `config.json` is shown.

***`We strongly recommend placing all repositories of tools / workflows in the same parent directory.`***

(All your repos should be side-by-side in sibling directories, as shown.)

```
......
"search_paths_cwl": {
        "global": [
            "../workflow-inference-compiler/cwl_adapters",
            "../image-workflows/cwl_adapters",
            "../biobb_adapters/biobb_adapters",
            "../mm-workflows/cwl_adapters"
        ],
        "gpu": [
            "../mm-workflows/gpu"
        ]
    },
"search_paths_yml": {
        "global": [
            "./workflow-inference-compiler/docs/tutorials",
            "../image-workflows/workflows",
            "../mm-workflows/examples"
        ]
    }
.....
```

If you do not specify config file using the command line argument `--config`, it will be automatically created for you the first time you run wic in `~/wic/global_config.json`. (Because of this, the first time you run wic you should be in the root directory of any one of your repos.) Then you can manually edit this file with additional sources of tools / workflows.

To avoid dealing with relative file paths in YAML files, by default

***`all tools / workflow names are required to be unique!`***

See [namespaces](advanced.md#namespaces) for details.

## Edges

What do the edges in a workflow represent? In many workflow languages (e.g. Argo), the edges represent dependencies between entire steps. Note that there could be ***`multiple`*** files or directories implicitly passed between two steps, but these workflow languages only model that as a single edge.

CWL models edges differently. In CWL, edges represent dependencies between ***`individual`*** explicit inputs and outputs. This fine-grained approach has several benefits, first and foremost increased parallelism. CWL also allows individual inputs and outputs to be tagged with metadata such as `type:` and `format:` tags. This additional information is what makes edge inference possible!

## Edge Inference Algorithm

First of all, a reminder that we can only connect an input in the current step to an output that already exists from some previous step.

The edge inference algorithm is actually rather simple: For each input with a given type and format, it checks for outputs that have the same type and format from one of the previous steps. Since many workflows are linear-ish pipelines, the steps are checked in reverse order (and the outputs of each step are also checked in reverse order). If there is a unique match, then great! If there are multiple matches, it arbitrarily chooses the first (i.e. most recent) match. For technical reasons edge inference is far from unique, so ***`users should always check that edge inference actually produces the intended DAG`***.

## Explicit Edges

If for some reason edge inference fails, you can always explicitly specify the edges using `'&var'` and `'*var'` notation. Simply use `'&var'` to create a reference to an output filename and then, in an input in any later step, use `'*var'` to dereference the filename and create an explicit edge between the output and the input. See [`examples/gromacs`](https://github.com/PolusAI/mm-workflows/blob/main/examples/gromacs) in `mm-workflows` repository for a concrete example. Due to yaml's [anchors and aliases](https://support.atlassian.com/bitbucket-cloud/docs/yaml-anchors/) notation (which you can still use!), these variables will need to be in quotes. (The notation is intended to be nearly identical, but instead of using `'*var'` to refer to the *contents* of `'&var'` it refers to the *path* to `'&var'`.)

## Inline Inputs

Of course, if you supply an input value directly, then the algorithm doesn't need to do either inference or explicit edges. The `message` input is a great example of this.

```
steps:
- echo:
    in:
      message: Hello World
```

Note that this is one key different between WIC and CWL. In CWL, all inputs must be given in a separate file. In WIC, inputs can be given inline and after compilation they will be automatically extracted into the separate file.