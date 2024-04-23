# Main Features / Design Overview

* Edge Inference

A large part of the complexity of workflows involves specifying the connections between steps. Every other workflow platform requires the user to explicitly specify the connections between steps, either graphically (e.g. KNIME) or via some verbose text format (e.g. CWL). This is fine for relatively trivial workflows, but rapidly becomes problematic for even mildly complicated workflows. Instead, we can automatically infer (almost) all of the edges! In addition to making life much easier on the user, it greatly helps to achieve the following high level goals.

* Subworkflows

It is critical that the steps of a workflow can be arbitrarily partitioned into logically related building blocks. Conversely, we must allow subworkflows to be arbitrarily composed together to make more complex workflows. Moreover, this must be done recursively, so that users can work at ever-higher levels of abstraction.

* Portable / Reproducible

Like other workflow platforms, a primary goal is to replace ad-hoc software tools (i.e. Bash scripts) with something that abstracts away machine-specific details (installed software packages, directory structure, etc) and can be easily shared with other researchers and reproduced.

There are several broad categories of workflow platforms; interested readers should check out the excellent introduction to the [SnakeMake paper](https://f1000research.com/articles/10-33/v2).

* Scalable

Like CWL, WIC workflows can be executed on a variety of platforms, from a laptop up to parallel clusters, cloud computing environments (e.g. AWS), and high performance computing environments.

* Declarative (No Programming Required!)

Our YAML format inherits all of the advantages of other declarative workflow platforms. Declarative formats are more restrictive than writing workflows in a programming language (e.g. Python), but this is actually an advantage because it facilitates additional tooling. For example:

* Graphical Workflow Builder

There is a work-in-progress effort to create a [WYSIWYG](https://en.wikipedia.org/wiki/WYSIWYG) GUI for building WIC workflows.

* No Walled Garden

Unlike some other workflow platforms which tend to lock-in users to their programming paradigm / APIs, we again follow the lead of the declarative platforms which allow arbitrary command line tools to be used. Open-source and closed-source, etc are all equally supported. If you can type it into a regular command line terminal, it can be supported.

## Advanced Features

See [advanced features](advanced.md) for details, but:

* Static Dispatch

WIC supports [ad hoc polymorphism](https://en.wikipedia.org/wiki/Ad_hoc_polymorphism) via [static dispatch](https://en.wikipedia.org/wiki/Static_dispatch). In other words, you can swap out 'interchangeable' implementations for any given step. Note that many other workflow platforms only allow swapping out implementations globally, for the workflow as a whole.

* Program Synthesis

In addition to Edge Inference, WIC has a (very limited) ability to partially write your workflows for you! See [program synthesis](advanced.md#program-synthesis)

* Realtime Monitoring

Users should be able to monitor the progress of an individual workflow step, even while it is still executing. This is important for steps which may take a long time to finish (e.g. supercomputer simulations). This is particularly useful for 'online' analysis algorithms.

To support this, WIC provides a way to invoke a CWL subinterpreter from any step within the main workflow. See [subinterpreters](advanced.md#subinterpreters)