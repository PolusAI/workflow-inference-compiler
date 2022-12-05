# Coding Standards

## Tooling

Tooling is essential to maximizing a developer's productivity and to minimizing the chance of introducing bugs. I highly recommend vscode, with the extensions listed in the [installation guide](../installguide.md#intellisense-code-completion). Most importantly, vscode provides [IntelliSense](https://code.visualstudio.com/docs/editor/intellisense) code completion.

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

### Code Coverage

We are using the pytest-cov code coverage plugin to make sure our tests are exercising most of the code. See `.coveragerc` for details. To generate a coverage report, simply use the `--cov` flag as shown below.

## CI/CD

Our Continuous Integration / Continuous Delivery files can be found in `.github/workflows/*.yml`. After every `git push`, this creates an isolated development environment and runs `mypy --no-incremental src/ tests/`, `pylint src/ tests/`, and `pytest --cov --workers 4`. Before pushing, please run `pytest --cov --workers 4 -m 'not slow'` or preferably the full `pytest --cov --workers 4`.

## Linting

We use pylint to check the code for style, formatting, and common mistakes. See `.pylintrc` for our configuration.

### Line Lengths

We are currently using `max-line-length=120`, although there are still some cases where lines exceed that length. Some developers feel strongly about a hard maximum of 100 or even 80 characters; we can try to accomodate a lower limit as necessary.
