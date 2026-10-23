"""A pytest plugin that makes plugin discovery an error.

Loaded with `-p core._poison_plugins`. Any module that reaches
`get_tools_cwl` or `get_yml_paths` under this plugin fails loudly instead of
quietly reading whatever the machine happens to have.
"""
from typing import Any, NoReturn

import sophios.plugins


def _poisoned(*_args: Any, **_kwargs: Any) -> NoReturn:
    raise AssertionError(
        'plugin discovery reached under the hermeticity run')


# setattr, not `sophios.plugins.get_tools_cwl = _poisoned`: pylint reads that assignment, takes every call to
# these functions in the repo for a call that never returns (NoReturn) and reports the code after it as unreachable.
for _name in ('get_tools_cwl', 'get_yml_paths'):
    setattr(sophios.plugins, _name, _poisoned)
