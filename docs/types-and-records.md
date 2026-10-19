# Types and records

Every input and output has a type, and Sophios uses types to catch a
mistaken connection before anything runs. In this chapter you will see what
Sophios checks, and what it leaves to the runner. Then you will meet the
**step-input record**, which you used for a condition in the last chapter,
and use it to feed one input from two steps at once.

## Types

A port's type is whatever its tool declares in CWL: a single kind of value
(`boolean`, `int`, `float`, `string`, `File`, `Directory` and a few more),
an array of one (`File[]`), a choice between several (`File?` is "a file or
nothing"), an enum, a record, or `Any`. Sophios has no types of its own; it
reads CWL's.

When a value meets a port, Sophios asks one question: can this value ever be
of the port's type?

- A **value you write** must have the port's type, with the three conversions
  [Inputs and values](inputs-and-values.md) describes. Anything else is
  `wic020`, in both forms.
- A **connection**, or a workflow input, is rejected only when Sophios can
  prove the two types never meet: two different kinds of value, an array
  against a single value, or choices with nothing in common. That is
  `wic023` in `.wic` and `api003` in Python, raised at the line that made
  the connection. When it cannot tell, because one side is `Any`, a record,
  an enum or a named type, it leaves the question to the runner, which
  checks the compiled workflow when it loads it.
- **Edge inference** is stricter than either, because it chooses rather than
  checks: an output is a candidate only when its type matches the input's
  ([Edge inference](edge-inference.md)).

So a `string` output connected to a `File` input is caught when you compile,
and an `Any` output connected to a `File` input is not caught by Sophios or
by cwltool's check of the compiled workflow; it fails only when the step
runs, if the value is not a file. Give tool outputs their exact type, and connect an
`Any` input by hand. A workflow input takes the type you declare under
`inputs:` in `.wic`, or the type of the first port it feeds in Python. A
workflow output takes the type of the output it names.

## Feed one input from two steps

Some tools take several files at once: a tool that joins files, merges
tables or compares samples. Its input is an array, `File[]`, and the files
come from several steps. Here two steps each make a file, and a third joins
them with `cat`:

```{literalinclude} cookbook/record.wic
:language: yaml
:caption: record.wic
```

```{literalinclude} cookbook/record.py
:language: python
:caption: record.py
```

The joining tool, `concat`, is written inline, under `run:` in `.wic` and as
a mapping given to `Step.from_cwl_document()` in Python; its one input,
`files`, takes a `File[]`. The input's value is the record: in `.wic`, `!cwl`
on a mapping; in Python, a `StepInput`. Its `source` lists the two outputs,
and `linkMerge: merge_flattened` joins them into one list. Run either form:

```bash
sophios --yaml record.wic --run_local
```

```bash
python record.py
```

Both succeed, and `joined.txt` in the output directory reads `Hello`: the
empty file followed by the line `echo` wrote. One file is empty on purpose,
so that the result shows plainly what came from where; with
`merge_flattened` the two files arrive as one list, in the order `source`
names them, and `cat` writes them out in that order.

## What a record can say

A record is CWL's step input written out in full. Its fields are CWL's, in
CWL's spelling in `.wic` and in Python's spelling in `StepInput`:

| CWL field | `StepInput` | Meaning |
|---|---|---|
| `source` | `source` | one output or workflow input, or a list of them |
| `linkMerge` | `link_merge` | how several sources combine: `merge_nested`, one list entry per source (the default), or `merge_flattened`, arrays joined into one |
| `pickValue` | `pick_value` | which sources to keep when some are `null`: `first_non_null`, `the_only_non_null` or `all_non_null` |
| `valueFrom` | `value_from` | an expression over the value, `self`, whose result is used instead |
| `default` | `default` | the value when the source gives none |
| `loadContents` | `load_contents` | read the first 64 KiB of a file into its `contents` |
| `loadListing` | `load_listing` | how deep a directory is listed |
| `label` | `label` | a title |

You need a record for five things: several sources into one input, as
above; a value for a `when:` under a name the tool does not take, as in
[Scatter and conditions](scatter-and-conditions.md); a `default`; a
transformation on the way in, such as
`minval: !cwl {source: count, valueFrom: '$(self + 1)'}`; and choosing among
the outputs of conditional steps with `pickValue`. For a single value or a
single connection, the plain spellings you already know are clearer.

A few rules keep records honest. A `source` holds connections only: a `!ii`
or a `!cwl` inside it, or a string in a `StepInput`'s `source`, is
`wic038` or `api001`. A record may bind a name the tool does not declare,
which is how a `when:` reads an extra value; that input's type is `Any`. A
record whose sources are merged, picked or transformed is not checked
against the input's type, since the transformation decides it. A key that is
not one of the fields above is `wic038`, and so is a record with no
`source`, `default` or `valueFrom`. [§3.5 of the `.wic`
reference](wic-reference.md#35-step-input-records-cwl-on-a-mapping) has the
rest.

## Recap

Sophios rejects a value whose type is wrong and a connection whose types can
never meet, and leaves to the runner the cases it cannot decide. A
step-input record, `!cwl {...}` in `.wic` or `StepInput(...)` in Python,
says more than one value can: several sources and how to merge them, a
default, a transformation, or a pick among conditional outputs.

## Next

Your pipeline has grown. In [Subworkflows](subworkflows.md) you will take
the part that makes the greeting and use it as a single step of a larger
workflow.
