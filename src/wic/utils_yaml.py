from typing import Any, Dict, Type

import yaml

# NOTE: In the following constructors, you CANNOT return the exact same yaml tag
# Otherwise, the loader is not idempotent. Specifically, then other tooling
# (i.e. the python api) cannot simply emit the dictionaries returned here,
# because then these constructors will fire again.


def anchor_constructor(loader: yaml.SafeLoader, node: yaml.nodes.ScalarNode) -> Dict[str, Dict[str, Any]]:
    key = loader.construct_scalar(node)
    name = 'wic_anchor'  # NOT '!&'
    return {name: {'key': key}}


def alias_constructor(loader: yaml.SafeLoader, node: yaml.nodes.ScalarNode) -> Dict[str, Dict[str, Any]]:
    key = loader.construct_scalar(node)
    name = 'wic_alias'  # NOT '!*'
    return {name: {'key': key}}


def inlineinput_constructor(loader: yaml.SafeLoader, node: yaml.nodes.ScalarNode) -> Dict[str, Dict[str, Any]]:
    val = loader.construct_scalar(node)
    name = 'wic_inline_input'  # NOT '!ii'
    return {name: {'val': val}}


def wic_loader() -> Type[yaml.SafeLoader]:
    loader = yaml.SafeLoader
    loader.add_constructor("!&", anchor_constructor)
    loader.add_constructor("!*", alias_constructor)
    loader.add_constructor("!ii", inlineinput_constructor)
    return loader
