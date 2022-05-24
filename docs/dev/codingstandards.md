# Coding Standards

## Tooling

Tooling is essential to maximizing a developer's productivity and to minimizing the chance of introducing bugs. I highly recommend vscode, with the extensions listed in the [installation guide](../installguide.md#vscode). Most importantly, vscode provides [IntelliSense](https://code.visualstudio.com/docs/editor/intellisense) code completion.

## MyPy Type Annotations

The code makes extensive use of mypy type annotations for static analysis. Together with the mypy vscode extension, this gives developers IntelliSense code completion which is absolutely critical. For various reasons mypy isn't perfect, but it's infinitely better than nothing at all. Please read the detailed comments in `wic_types.py` for additional considerations w.r.t. [Algebraic Data Types](https://en.wikipedia.org/wiki/Algebraic_data_type).

## Docstrings

I use the vscode extension `autoDocstring` to help generate docstrings. If you wait until after you add mypy type annotations to all of a function's arguments, autoDocstring will insert the types into the docstrings as well. Docstrings complement mypy type annotations; they do not replace them. We need to use both.

The contents of the docstrings should try to describe the purpose of the function and the arguments at a high level. For example:

```python
def step_name_str(yaml_stem: str, i: int, step_key: str) -> str:
    """Returns a string which uniquely and hierarchically identifies a step in a workflow

    Args:
        yaml_stem (str): The name of the workflow (filepath stem)
        i (int): The (zero-based) step number
        step_key (str): The name of the step (used as a dict key)

    Returns:
        str: The parameters (and the word 'step') joined together with double underscores
    """
    ...
```

## Integration / Regression Tests

We have a small number of integration / regression tests located in `tests/`. These tests are very powerful and have been extremely helpful in catching various issues. The command `pytest` will run the full set of tests, which takes about 20-30 minutes on a laptop. For interactive development purposes this is too slow, so you can use `pytest -m fast` or `pytest -m 'not slow'`. These commands test the compiler only (not the runtime) and only take about 15 seconds and 2 minutes, respectively.

In addition to the above tests, we need to add some unit tests. Due to the highly recursive nature of the compilation algorithm, it has proven difficult to test functions in isolation. There will probably need to be some refactoring to facilitate this.

## CI/CD

Our Continuous Integration / Continuous Delivery files can be found in `.github/workflows/*.yml`. After every `git push`, this creates an isolated development environment and runs `pytest`. Before pushing, please run `pytest -m 'not slow'` or preferably the full `pytest`.

## Line Lengths

I personally don't mind excessively long lines at all, but others may feel more strongly on this. Side-by-side diffs certainly look better with a maximum line length, but the occasional line wrap doesn't bother me. That said, I use `"python.formatting.blackArgs": ["--line-length", "120"],` in my `.vscode/settings.json`
