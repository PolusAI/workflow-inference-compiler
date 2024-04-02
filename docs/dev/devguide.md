# Developer Guide

See [algorithms](algorithms.md) for a description of the compilation algorithms and some high-level implementation considerations. I hope you like recursion! ;)

## Coding Standards

See [coding standards](codingstandards.md)

## Git Etiquette

See [git etiquette](gitetiquette.md)

## Known Issues

### Bad User Inputs

Although we now have a formal YAML schema and perform validation, I'm sure there are plenty of other ways in which users can crash the system, so we need to make an effort to find these cases and do more error checking.

### Uniqueness and Dict keys

There are some cases where we blindly attempt to index into a dict with a suspicious key. If the key is not found, it will generate a nasty stack trace. We need to find all instances of not-so-good keys and do additional error checking.

### Globbing Unexpected File Order

If a workflow step generates a file with some association between a filename and numerical value, if the numerical values are being extracted in order of each row/column from the file then there is no guarantee that using glob to retrieve files will produce the same order consistent with the extracted numerical value array. Another example is illustrated below without scattering, where using glob will result in inconsistency between input order and output order of files. With scattering it is possible in some cases to induce the correct output ordering consistent with input order, however it is best practice to adopt reading filename indices from an output file rather than using glob to ensure consistent order. This way developer does not have to think about which cases glob might or might now work in.

```
input: [3, 2, 1]  # Here is the order of the input array.
```

```
#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
baseCommand: ["python", "example.py"]

requirements:
- class: InitialWorkDirRequirement
  listing:
  # See https://www.commonwl.org/user_guide/topics/creating-files-at-runtime.html
    - entryname: example.py
      entry: |
        import sys
        from pathlib import Path

        for arg in sys.argv[1:]:
            Path(f'{arg}.txt').touch()

inputs:
  input:
    type: Any[]
    inputBinding:
      position: 0

outputs:
  output:
    type: File[]
    outputBinding:
      glob: "*.txt"
```

```
cwltool touch_array.cwl touch_array_inputs.yml


INFO [job test.cwl] /tmp/0x0q86bg$ python \
    example.py \
    3 \
    2 \
    1
INFO [job test.cwl] completed success
{
    "output": [
        {
            "location": "file:///home/walkerbd/workflow-inference-compiler/cwl_adapters/1.txt",
            "basename": "1.txt",
            "class": "File",
            "checksum": "sha1$da39a3ee5e6b4b0d3255bfef95601890afd80709",
            "size": 0,
            "path": "/home/walkerbd/workflow-inference-compiler/cwl_adapters/1.txt"
        },
        {
            "location": "file:///home/walkerbd/workflow-inference-compiler/cwl_adapters/2.txt",
            "basename": "2.txt",
            "class": "File",
            "checksum": "sha1$da39a3ee5e6b4b0d3255bfef95601890afd80709",
            "size": 0,
            "path": "/home/walkerbd/workflow-inference-compiler/cwl_adapters/2.txt"
        },
        {
            "location": "file:///home/walkerbd/workflow-inference-compiler/cwl_adapters/3.txt",
            "basename": "3.txt",
            "class": "File",
            "checksum": "sha1$da39a3ee5e6b4b0d3255bfef95601890afd80709",
            "size": 0,
            "path": "/home/walkerbd/workflow-inference-compiler/cwl_adapters/3.txt"
        }
    ]
```

As can be seen from the output json blob the order returned is not the same as input order.