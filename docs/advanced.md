# Advanced Features

## Edge Inference Configuration

### Naming Conventions
If `--inference_use_naming_conventions` is enabled, matches can be refined based on the naming conventions of the inputs and outputs in the CWL CommandLineTools. Specifically, first `input_` is removed from the input name and `output_` is removed from all output names. Then, the default renamings contained in `renaming_conventions` tag of `config.json` (shown below) are iteratively applied to the input name (only), and then the modified input name and all of the output names are checked for equality.

```
"renaming_conventions": [
        [
            "energy_",
            "edr_"
        ],
        [
            "structure_",
            "tpr_"
        ],
        [
            "traj_",
            "trr_"
        ]
    ]
```

If there is now a unique match, then great! If there are still multiple matches, it chooses the first (i.e. most recent) match. If there are now no matches, it ignores the naming conventions and chooses the first (i.e. most recent) match based on types and formats only. If there are still multiple matches, it again chooses the first (i.e. most recent) match. Note that there are cases (i.e. file format conversions) where using naming conventions may not yield the desired behavior, so again ***`users should always check that edge inference actually produces the intended DAG`***.

## Static dispatch

WIC supports [ad hoc polymorphism](https://en.wikipedia.org/wiki/Ad_hoc_polymorphism) via [static dispatch](https://en.wikipedia.org/wiki/Static_dispatch).

To use static dispatch, first you need to create a YAML file which aggregates the implementation-specific workflows:

```yaml
wic:
  default_implementation: implementation1
  implementations:
    implementation1:
      steps:
        - implementation1.wic:
    implementation2:
      steps:
        - implementation2.wic:
```

Then you just need to choose a specific implementation at the call site:

```yaml
steps:
  ... some workflow steps
  - static_dispatch.wic:
  ... some more workflow steps

wic:
  steps:
    (2, static_dispatch.wic):
      wic:
        implementation: implementation2
```

The most common use case of static dispatch is to swap out 'identical' subworkflows. However, this constraint is intentionally not enforced and it is completely up to the user. In fact, you may want to swap out two implementations that use different algorithms to achieve the same high-level goal.

### Program Synthesis

WIC supports a very limited form of [Program Synthesis](https://en.wikipedia.org/wiki/Program_synthesis). The `--insert_steps_automatically` flag augments the inference algorithm by attempting to insert an extra step if inference initially fails. For now, the possible steps are limited to a whitelist of known-good steps whose names start with `insert_steps_automatically_`. For details on the algorithm, see [speculative compilation](dev/algorithms.md#speculative-compilation). In short, if inference fails, we can attempt to insert a step and try again.

NOTE: This should be considered an experimental feature; it is very useful for some specific use cases (e.g. file format conversions), but more work is needed to see how useful this particular algorithm is in general.

#### Known issues

Note however that while insert_steps_automatically_*.cwl files can come from pre-compiled wic subworkflows, it is currently necessary to find & replace all instances of triple underscores ___ with double underscores __. (This is because triple underscores are reserved/interpreted by the compiler as ‘internal’ namespaceing, and in this case we want to treat pre-compiled wic files as a black box. See the [dev guide](dev/algorithms.md#namespacing) for the gory details.)

## Subinterpreters

For realtime monitoring, we want to asynchronously run an auxiliary workflow while the main workflow is still running. We implement this by invoking a subinterpreter.

The `cwl_subinterpreter` subinterpreter will repeatedly run an auxiliary workflow which is completely independent from the main workflow, upto some fixed number of iterations. (A fixed number of iterations is used because the main interpreter and subinterpreter are completely independent; there is no way of passing status information between the interpreters.)

NOTE: This should be considered an experimental feature; the CWL standard does not consider realtime monitoring (and/or other implementation-specific details).

## YAML Metadata Annotations

For various reasons, one may want to extend a CWL file with additional tags. We have decided to jam all of these additional tags into a single top-level `wic:` tag, so that only one non-CWL tag needs to be removed during compilation.

YAML files can be annotated with metadata inside of a top-level `wic:` tag. Metadata that applies to the entire workflow can be specified within the first `wic:` tag. To apply metadata to individual steps, use a `steps:` tag and then specify the step using (step_num, tool_name) where step_num is 1-based. (Since a given tool can be used more than once in a workflow, a step number is necessary to uniquely identify the desired step.) Then, add the desired metadata under another `wic:` tag.

### GraphViz options

A simple use case is providing labeling and alignment metadata for generating the Graphviz DAGs. There is also a `style: invis` tag, which can be used to hide certain nodes. For now, these three tags are sufficient for generating visually appealing DAGs. We do not anticipate needing to support many additional graphviz features.

```yaml
...
wic:
  graphviz:
    label: Descriptive Subworkflow Name
    ranksame:
    - (1, short_step_name_1)
    - (5, short_step_name_5)
  steps:
    (1, short_step_name_1):
      wic:
        graphviz:
          label: Descriptive Step Name 1
...
```

### Overloading / Parameter Passing

Since every `*.wic` file may contain a `wic:` tag, we can recursively pass in parameters / recursively overload metadata. Values from parent `*.wic` files overwrite values in the child `*.wic` files. Note that we do not need to modify any of the child workflows!

Also note that this is done statically, at compile time (not runtime).

Thus we retain edge inference, explicit edges, composability, reusability, and have even gained customizability!

```yaml
...
# Put everything under one top-level wic: tag to facilitate easy merging and removal.
wic:
  steps:
    (1, some_subworkflow.wic):
      wic:
        steps:
          (2, another_subworkflow.wic):
            wic:
              steps:
                (3, step_name):
                  in:
                    input_name: new_input_value
                  out:
                  - output_name: !& new_anchor_name
...
```

## Namespaces

Namespaces can be used to distinguish two different tools / workflows with the same name from different sources. For example, suppose a collaborator has shared an alternative minimization protocol, which we have downloaded to `bar/min.wic`. We can use their protocol by adding the namespace tag `foo` to `search_paths_wic` tag of `config.json` and annotating the call site with `namespace: foo` as shown below.

```yaml
...
wic:
  steps:
    (1, min.wic):
      wic:
        namespace: foo
...
```

### Known Issues

Note that within a single namespace, names are still required to be unique. If your CWL and/or YAML filenames are *not* unique, the message `Warning: overwriting <filestem> in namespace <namespace>` will be printed during compilation, referring to the fact that the in-memory representation (only) is being overwritten. This will likely cause severe problems, so please heed this particular warning!

On the other hand, that message may also be printed when performing [speculative compilation](dev/algorithms.md#speculative-compilation), but in that case it is expected and should not be considered a warning.

## Miscellaneous

### Inference Rules

This section really belongs in the Edge Inference Configuration section, but it is less important.

Users can customize the inference algorithm using inference rules. The default inference rules stored in `inference_rules` tag of `config.json` are shown below:

```
"inference_rules": {
        "edam:format_3881": "continue",
        "edam:format_3987": "continue",
        "edam:format_3878": "break",
        "edam:format_2033": "break"
    }
```

Currently, the only inference rule implemented is `break`, which stops the inference algorithm from considering any further outputs beyond the current output from matching the current input. (The current output is allowed, i.e. break is inclusive.) This is useful when the most recent output file is desired, but the inference algorithm for some reason doesn't match it and chooses a subsequent / earlier file. This can happen when converting from one file format, performing a workflow step, and converting back to the original format, where in some cases the inference algorithm may choose the original file, thus accidentally skipping the workflow step.

The next rule to be added will be `continue`, which will simply skip over the current output and continue processing the remaining outputs. This is useful when the 'original' version of a file should be preferred. Again, when performing file format conversions, some auxiliary files (i.e. files which are needed in the intermediate workflow step but are not otherwise intended to be modified) may be accidentally modified due to lossy file format conversions. If we know that a double-conversion ought to be the identity but we do not want to rely on that, we can use `continue` to invalidate the double-converted file and thus the inference algorithm should eventually match on the original file.