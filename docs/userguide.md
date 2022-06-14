# User Guide

Recall for a moment the vague instructions your PhD advisor hastily scribbled onto the chalkboard about how to do a calculation. Now imagine that those scribbles were actually executable! That's the goal! The goal is to allow high-level, domain-specific concepts to be directly specified in a user-friendly YAML format. Then the abstract scientific protocol is automatically translated into specific concrete steps, executed on a remote job cluster, and automated analysis is performed.

## Main Features / Design Overview

See [overview](overview.md)

## Plugins

Before we get started, we will need to configure the compiler with plugins. ('Plugins' are either CWL CommandLineTools or YAML workflow files.) By default, the compiler will search for plugins within the directories and subdirectories listed in `cwl_dirs.txt` and `yml_dirs.txt`. To enable additional plugins, simply download the plugins and add additional lines to these files. 

```
# Namespace Directory
global      biobb/
global      cwl_adapters/
# foo       a/relative/path/
# bar       /an/absolute/path/
```

To avoid having to deal with relative file paths in YAML files, all plugin names (e.g. `mdrun` & `gmx_energy` below) are required to be unique. See [plugin namespaces](userguide.md#plugin-namespaces) for details.

## Edge Inference Algorithm

The edge inference algorithm is actually rather simple: For each input in the current step of a workflow, it checks for compatible outputs in the previous steps. Since most operations presumably desire the most recent compatible output, by default the steps are checked in reverse order and the outputs of each step are also checked in reverse order. For each output, it first checks for matching types and formats. If there is a unique match, then great! If there are multiple matches, it tentatively chooses the first (i.e. most recent) match. For technical reasons edge inference is far from unique, so ***`users should always check that edge inference actually produces the intended DAG`***.

### Naming Conventions
If `--cwl_inference_use_naming_conventions` is enabled, matches can be refined based on the naming conventions of the inputs and outputs in the CWL CommandLineTools. Specifically, first `input_` is removed from the input name and `output_` is removed from all output names. Then, the default renamings contained in `renaming_conventions.txt` (shown below) are iteratively applied to the input name (only), and then the modified input name and all of the output names are checked for equality.

```
# The biobb CWL files do not always use consistent naming
# conventions, so we need to perform some renamings here.
# Eventually, the CWL files themselves should be fixed.

energy_     edr_
structure_  tpr_
traj_       trr_
```

If there is now a unique match, then great! If there are still multiple matches, it chooses the first (i.e. most recent) match. If there are now no matches, it ignores the naming conventions and chooses the first (i.e. most recent) match based on types and formats only. If there are still multiple matches, it again chooses the first (i.e. most recent) match. Note that there are cases (i.e. file format conversions) where using naming conventions may not yield the desired behavior, so again ***`users should always check that edge inference actually produces the intended DAG`***.

### Inference Rules

Users can customize the inference algorithm using inference rules. The default inference rules stored in inference_rules.txt are shown below:

```
# Amber, gromacs (zipped 3880) topology
edam:format_3881 continue
edam:format_3987 continue

# Amber, gromacs coordinates
edam:format_3878 break
edam:format_2033 break
```

Currently, the only inference rule implemented is `break`, which stops the inference algorithm from considering any further outputs beyond the current output from matching the current input. (The current output is allowed, i.e. break is inclusive.) This is useful when the most recent output file is desired, but the inference algorithm for some reason doesn't match it and chooses a subsequent / earlier file. This can happen when converting from one file format, performing a workflow step, and converting back to the original format, where in some cases the inference algorithm may choose the original file, thus accidentally skipping the workflow step.

The next rule to be added will be `continue`, which will simply skip over the current output and continue processing the remaining outputs. This is useful when the 'original' version of a file should be preferred. Again, when performing file format conversions, some auxillary files (i.e. files which are needed in the intermediate workflow step but are not otherwise intended to be modified) may be accidentally modifed due to lossy file format conversions. If we know that a double-conversion ought to be the identity but we do not want to rely on that, we can use `continue` to invalidate the double-converted file and thus the inference algorithm should eventually match on the original file. The prototypical example of such a file is a molecular 'topology' file; the coordinates will obviously change after each molecular dynamics simulation, but the topology represents hardcoded connectivity information which should not change.

## Explicit Edges

If for some reason edge inference fails, you can always explicitly specify the edges using `'&var'` and `'*var'` notation. Simply use `'&var'` to create a reference to an output filename and then, in an input in any later step, use `'*var'` to dereference the filename and create an explicit edge between the output and the input. See examples/gromacs for a concrete example. Due to yaml's [anchors and aliases](https://support.atlassian.com/bitbucket-cloud/docs/yaml-anchors/) notation (which you can still use!), these variables will need to be in quotes. (The notation is intended to be nearly identical, but instead of using `'*var'` to refer to the *contents* of `'&var'` it refers to the *path* to `'&var'`.)

## Backend Independence

To use backend independence, first you need to create a YAML file which aggregates the backend-specific workflows:

```yaml
wic:
  default_backend: gromacs
  backends:
    gromacs:
      steps:
        - npt_gromacs.yml:
    amber:
      steps:
        - npt_amber.yml:
  graphviz:
    label: Constant Pressure
```

Then you just need to choose a specific backend at the call site:

```yaml
steps:
  - nvt.yml:
  - npt.yml:

wic:
  graphviz:
    label: Equilibration
  steps:
    (2, npt.yml):
      wic:
        backend: amber
```
This will override the default backend of `gromacs` and use `amber`. This really just means that `npt_amber.yml` is called instead of `npt_gromacs.yml` The system will automatically insert the necessary file format conversions as determined below.

The most common use case of backend independence is to swap out 'identical' subworkflows. However, this constraint is intentionally not enforced and it is completely up to the user. In fact, you may want to swap out two backends that use different algorithms to achieve the same high-level goal. In other words, 'backend independence' is really just a form of [ad hoc polymorphism](https://en.wikipedia.org/wiki/Ad_hoc_polymorphism) and/or [dynamic dispatch](https://en.wikipedia.org/wiki/Dynamic_dispatch).

### Automated File Format Conversion / Speculative Compilation

When using backend independence, it will often require file format conversion steps. The compiler can now automatically insert these steps into the workflow! For details on the algorithm, see [speculative compilation](dev/algorithms.md#speculative-compilation). In short, if inference fails, we insert a file format conversion and try again.

For now, the possible file format conversions are limited to a whitelist of known good steps whose names start with `conversion_`. However, there is no fundamental reason to limit this speculative complation strategy to only file format conversions. In fact, we can automatically insert arbitrary subworkflows!

#### Known issues

Note however that while conversion_*.cwl files can come from pre-compiled yml subworkflows, it is currently necessary to find & replace all instances of triple underscores ___ with double underscores __. (This is because triple underscores are reserved/interpreted by the compiler as ‘internal’ namespaceing, and in this case we want to treat pre-compiled yml files as a black box. See the [dev guide](dev/algorithms.md#namespacing) for the gory details.)

## Real-time analysis / Speculative Execution

Ordinarily the runtime system will wait until the previous step(s) are all complete before executing the next step. However, for real-time analysis we want to speculatively execute an arbitrary subworkflow (i.e. to parse log files, etc) before the the previous step(s) have finished. (Note that 'previous' is w.r.t. the DAG, i.e. it refers to all of the nodes which are dependencies of the current step.)

Speculative execution is currently implemented by `cwl_watcher`, which invokes a second instance of the runtime system separately and asynchronously. A portion of `examples/gromacs/nvt.yml` is shown below. You can see that the `in:` tag of gmx_energy is identical to the `config:` tag of cwl_watcher. This currently needs to be manually copy & pasted (and indented), but it should be possible to automatically do this in the future.

```yaml
...
  - mdrun:
      in:
        output_edr_path: '&nvt.edr' # Explicit edge reference / anchor
        # (This edge can be inferred, but made explicit for demonstration purposes.)
  - gmx_energy:
      in:
        input_energy_path: '*nvt.edr' # Explicit edge dereference / alias
        config:
          terms: [Temperature]
        output_xvg_path: temperature.xvg
# NOTE: explicit edges are not supported with cwl_watcher, and all filenames
# must be globally unique!
  - cwl_watcher:
      in:
        #cachedir_path: /absolute/path/to/cachedir/ (automatically filled in by wic)
        file_pattern: '*nvt.edr' # This * is a glob wildcard, NOT an explicit edge!
        cwl_tool: gmx_energy # This can also be an arbitrary subworkflow!
        max_times: '5'
        config:
          in:
            input_energy_path: '*nvt.edr' # This * is automatically removed.
            config:
              terms: [Temperature]
            output_xvg_path: temperature.xvg
...
```

Note that although gmx_energy appears before cwl_watcher in the YAML file, gmx_energy is independent of cwl_watcher in the DAG and thus not considered to be a previous step. We include gmx_energy simply to guarantee that the analysis gets run one more time in the main workflow, when all the files are known to be in their final state.

### Known Issues

Since the two runtimes are not linked, there is not currently a reliable way to determine if the previous steps have finished. Thus, to guarantee termination of the second runtime, we simply execute `cwl_tool` upto `max_times`. We also waive any guarantees about the files, so the subworkflow in the second runtime may of course fail for any number of reasons. Thus, we do not propagate speculative failures up to the main workflow.

The runtime system intentionally hides the working sub-directories of each step. Thus, we are forced to use a file watcher (hence the name cwl_watcher) recursively starting from `cachedir_path`. This is why all filenames used with cwl_watcher must be globally unique. (Actually, for technical reasons we cannot use a file watching library; we simply use a good old fashioned polling loop.)

## Real-time plots

It is assumed that the real-time analysis takes care of the complex log file parsing, etc and produces simple tabular data files (i.e. csv files separated by whitespace instead of a comma). We need to use the same file watching / polling trick as above to locate these tabular data files. The first argument to the following command is the directory in which to look for the files. (By default it is `cachedir` because that is the default value of the  `--cachedir` wic command line argument.) You can also optionally supply the file patterns, which by default are `*.xvg` and `*.dat`.

```
python RealtimePlots.py cachedir
```

## Labshare Compute

As previously mentioned, one of the beautiful things about the declarative approach to workflows is that we can execute workflows on massive machines just as easily as executing workflows on a local laptop. Concretely, merely changing `--cwl_run_local` to `--cwl_run_slurm`, we can execute the exact same workflow on the NCATS HPC cluster! That's it! Absolutely no modifications necessary!

### Authentication Access Token

When using `--cwl_run_slurm` you will also need to use `--compute_access_token $ACCESS_TOKEN`. Unfortunately, there is currently no programmatic way of obtaining the access token via an API. You will need to manually perform the following steps:

* Go to https://a-qa.labshare.org
* Click Login and then Azure
* Authenticate (using Microsoft Authenticator)
* You should be redirected to https://a-qa.labshare.org/tenants/polus-qa/dashboard
* Open the Javascript Console
* Click on Session Storage
* Click on the key oidc.user:https://a-qa.labshare.org/_auth/vTest/auth/polus-qa/:auth-ui
* Copy the massive hash string under the "access_token" tag
* In a bash terminal, create the environment variable `export ACCESS_TOKEN=...` (where ... means paste in the hash string)

Unfortunately, the access_token currently expires after about an hour (see "expires_at" tag), so you will need to repeat these steps periodically.

### Workflow Status

After submitting a workflow, users can check on its status by logging into login.hpc.ncats.io and using the `squeue` command. Currently, the workflows are executed under the placeholder svc-polus user.

```
ssh login.hpc.ncats.io
```

```
watch -n5 squeue -u svc-polus
```

After a ~2-3 minute initial delay (due to provisioning, etc) you should see nodes starting, running, and finishing, corresponding to the individual steps in the workflow. The output files and logs are currently stored under /project/labshare-compute/

## YAML Metadata Annotations

YAML files can be annotated with metadata inside of a top-level `wic:` tag. Metadata that applies to the entire workflow can be specified within the first `wic:` tag. To apply metadata to individual steps, use a `steps:` tag and then specify the step using (step_num, tool_name) where step_num is 1-based. (Since a given tool can be used more than once in a workflow, a step number is necessary to uniquely identify the desired step.) Then, add the desired metadata under another `wic:` tag.

### GraphViz options

A simple use case is providing labeling and alignment metadata for generating the Graphviz DAGs. There is also a `style: invis` tag, which can be used to hide certain nodes. For now, these three tags are sufficient for generating visually appealing DAGs. We do not anticipate needing to support many additional graphviz features.

A portion of `examples/gromacs/setup.yml` is shown below.

```yaml
...
wic:
  graphviz:
    label: System Setup
    ranksame:
    - (1, pdb2gmx)
    - (5, genion)
  steps:
    (1, pdb2gmx):
      wic:
        graphviz:
          label: 'Generate\nInitial\nTopology'
...
```

### Overloading / Parameter Passing

This example shows how we can recursively pass in parameters / recursively overload metadata.

Suppose we want to do a very careful minimization, first in vacuum and then in solvent (i.e. `examples/gromacs/setup_vac_min.yml`). We would like to re-use the abstract minimization protocol from `min.yml`. However, our stability analysis requires an explicit edge definition from the final minimized coordinates (i.e. in solvent). If we try to simply add `output_tpr_path: '&min.tpr'` directly to `min.yml`, there will be duplicate definitions! This is not allowed (it will generate an exception).

The solution is to pass in this parameter to only the second instance of `min.yml`. Since every `*.yml` file may contain a `wic:` tag, this is implemented by simply recursively merging the dictionaries, where values from parent `*.yml` files overwrite values in the child `*.yml` files. Note that we do not need to modify `min.yml`!

Thus we retain edge inference, explicit edges, composability, reusability, and have even gained customizability!

A portion of `examples/gromacs/basic.yml` is shown below.

```yaml
...
# Put everything under one top-level wic: tag to facilitate easy merging and removal.
wic:
  graphviz:
    label: Molecular Dynamics
  steps:
    (1, min.yml):
      wic:
        steps:
          (2, cg.yml):
            wic:
              steps:
                (1, grompp):
                  in:
                    output_tpr_path: '&min.tpr'
...
```

## Plugin Namespaces

Namespaces can be used to distinguish two different plugins with the same name from different sources. For example, suppose a collaborator has shared an alternative minimization protocol, which we have downloaded to `bar/min.yml`. We can use their protocol by adding the line `foo    bar/` to `yml_dirs.txt` and annotating the call site in `basic.yml` with `namespace: foo` as shown below.

```yaml
...
wic:
  graphviz:
    label: Molecular Dynamics
  steps:
    (1, min.yml):
      wic:
        namespace: foo
...
```

### Known Issues

Note that within a single namespace, names are still required to be unique. If your CWL and/or YAML filenames are *not* unique, the message `Warning: overwriting <filestem> in namespace <namespace>` will be printed during compilation, referring to the fact that the in-memory representation (only) is being overwritten. This will likely cause severe problems, so please heed this particular warning!

On the other hand, that message may also be printed when performing [speculative compilation](dev/algorithms.md#speculative-compilation), but in that case it is expected and should not be considered a warning.