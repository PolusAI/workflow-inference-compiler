# git etiquette

## Branching Model

I am generally in favor of the [Scaled Trunk-Based Development](https://trunkbaseddevelopment.com) model. In this model, branches should be short-lived (a few days to a week) and frequently merged back into master/main after all tests pass. This model tends to avoid massive merge conflicts, which are no fun. This is particularly important in the early stages of a project when there are major refactorings happening and commits do not always commute.

## rebaseing 'vs' squashing

Ideally, branches should be cleaned up before merging using `git rebase --interactive`. Rebaseing can be incredibly powerful, but rewriting history can also be incredibly dangerous. Thus, since many users rarely rebase before merging, I am generally in favor of using `git merge --squash`. The history is still available in the other branch, but this gives master/main a cleaner history.

## Commit Frequency

In general, smaller commits are better, at least from the point of view that is easier to squash small commits together than [splitting commits](https://git-scm.com/docs/git-rebase#_splitting_commits) into individual hunks. Commits should try to be self-contained, and only change one logically related feature at a time.

### Pull Requests

Like commits, PRs should try to be self-contained and small and most importantly easily reviewable. If your PR is greater than ~250 lines, the probability of getting it merged drops rapidly.

## Testing and CI/CD

Before pushing, please run `pytest -m serial` or preferably the full `pytest -m serial && pytest -m "not serial" --workers 4`.

## .gitignore

We have a `.gitignore` file to prevent developers from accidentally pushing certain files. However, .gitignore files are somewhat fragile and should not really be relied upon. In other words, please avoid `git add *` and please use `git status` to check that you are *only* adding the specific files you intend to add.