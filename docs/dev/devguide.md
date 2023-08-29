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