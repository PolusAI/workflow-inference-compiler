# Developer Guide

See [algorithms](algorithms.md) for a description of the compilation algorithms and some high-level implementation considerations. I hope you like recursion! ;)

## Coding Standards

See [coding standards](codingstandards.md)

## Git Etiquette

See [git etiquette](gitetiquette.md)

## Git Submodules

The plugins (i.e. CWL adapters) are added to this repo as git submodules. This is because they should be in principle completely independent of the core inference and compilation algorithms and because they should not need to be modified very often.

As shown in the README, cloning the main repo will not clone the submodules; you will need to run the following command:
```
git submodule init && git submodule update
```

Developers should be very careful when using git submodules! The following links explain why:

* [how do i commit changes in a git submodule](https://stackoverflow.com/questions/5542910/how-do-i-commit-changes-in-a-git-submodule)
* [how to link git repos](https://stackoverflow.com/questions/36554810/how-to-link-folder-from-a-git-repo-to-another-repo)
* [git submodule core concepts](https://www.atlassian.com/git/articles/core-concept-workflows-and-tips)

"Why not use `git subtree`?" I'm not opposed to alternative git workflows, but I think the independent and static nature of the plugins is well-suited to submodules.

## Known Issues

### Bad User Inputs

Although we now have a formal YAML schema and perform validation, I'm sure there are plenty of other ways in which users can crash the system, so we need to make an effort to find these cases and do more error checking.

### Uniqueness and Dict keys

There are some cases where we blindly attempt to index into a dict with a suspicious key. If the key is not found, it will generate a nasty stack trace. We need to find all instances of not-so-good keys and do additional error checking.