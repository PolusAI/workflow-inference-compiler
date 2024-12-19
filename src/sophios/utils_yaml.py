from typing import Any, Dict, Type

import yaml

# NOTE: In the following constructors, you CANNOT return the exact same yaml tag
# Otherwise, the loader is not idempotent. Specifically, then other tooling
# (i.e. the python api) cannot simply emit the dictionaries returned here,
# because then these constructors will fire again.


def anchor_constructor(loader: yaml.SafeLoader, node: yaml.nodes.ScalarNode) -> Dict[str, Any]:
    val = loader.construct_scalar(node)
    name = 'wic_anchor'  # NOT '!&'
    return {name: val}


def alias_constructor(loader: yaml.SafeLoader, node: yaml.nodes.ScalarNode) -> Dict[str, Any]:
    val = loader.construct_scalar(node)
    name = 'wic_alias'  # NOT '!*'
    return {name: val}


def inlineinput_constructor(loader: yaml.SafeLoader, node: yaml.nodes.Node) -> Dict[str, Dict[str, Any]]:
    val: Any
    if isinstance(node, yaml.nodes.ScalarNode):
        try:
            # loader.construct_scalar always returns a string, whereas
            if node.value == "":
                val = ""
            else:
                val = yaml.safe_load(node.value)
            # yaml.safe_load returns the correct primitive types
        except Exception:
            # but fallback to a string if it is not actually a primitive type.
            val = loader.construct_scalar(node)
    elif isinstance(node, yaml.nodes.MappingNode):
        val = loader.construct_mapping(node)
    elif isinstance(node, yaml.nodes.SequenceNode):
        val = loader.construct_sequence(node)
    else:
        raise Exception(f'Unknown yaml node type! {node}')
    name = 'wic_inline_input'  # NOT '!ii'
    return {name: val}


def wic_loader() -> Type[yaml.SafeLoader]:
    loader = yaml.SafeLoader
    loader.add_constructor("!&", anchor_constructor)
    loader.add_constructor("!*", alias_constructor)
    loader.add_constructor("!ii", inlineinput_constructor)
    return loader
