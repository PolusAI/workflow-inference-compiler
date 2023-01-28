from typing import Any, Dict, List, Tuple, Type

import yaml

# snakeyaml (a cromwell dependency) refuses to parse yaml files with more than
# 50 anchors/aliases to prevent Billion Laughs attacks.
# See https://en.wikipedia.org/wiki/Billion_laughs_attack
# Solution: Inline the contents of the aliases into the anchors.
# See https://ttl255.com/yaml-anchors-and-aliases-and-how-to-disable-them/#override


class NoAliasDumper(yaml.SafeDumper):
    def ignore_aliases(self, data: Any) -> bool:
        return True


def anchor_constructor(loader: yaml.SafeLoader, node: yaml.nodes.MappingNode) -> Dict[str, Dict[str, Any]]:
    pairs: List[Tuple[str, Any]] = loader.construct_pairs(node)
    if not len(pairs) == 1:
        raise Exception("Error! " + str(dict(pairs)) + " must be a single {key, val}")
    (key, val) = pairs[0]
    return {'!&': {'key': key, 'val': val}}


def alias_constructor(loader: yaml.SafeLoader, node: yaml.nodes.ScalarNode) -> Dict[str, Dict[str, Any]]:
    key = loader.construct_scalar(node)
    return {'!*': {'key': key}}


def inlineinput_constructor(loader: yaml.SafeLoader, node: yaml.nodes.ScalarNode) -> Dict[str, Dict[str, Any]]:
    val = loader.construct_scalar(node)
    return {'!ii': {'val': val}}


def wic_loader() -> Type[yaml.SafeLoader]:
    loader = yaml.SafeLoader
    loader.add_constructor("!&", anchor_constructor)
    loader.add_constructor("!*", alias_constructor)
    loader.add_constructor("!ii", inlineinput_constructor)
    return loader
