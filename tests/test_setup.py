import argparse
import sys
import time
from typing import Dict
from unittest.mock import patch

import hypothesis_jsonschema as hj

import wic
import wic.cli
import wic.main
import wic.schemas
import wic.schemas.wic_schema
import wic.utils
from wic.wic_types import Json, Yaml


def get_args(yaml_path: str = '') -> argparse.Namespace:
    """This is used to get mock command line arguments.

    Returns:
        argparse.Namespace: The mocked command line arguments
    """
    testargs = ['wic', '--yaml', yaml_path, '--cwl_output_intermediate_files', 'True']  # ignore --yaml
    # For now, we need to enable --cwl_output_intermediate_files. See comment in compiler.py
    with patch.object(sys, 'argv', testargs):
        args: argparse.Namespace = wic.cli.parser.parse_args()
    return args


tools_cwl = wic.main.get_tools_cwl(get_args().cwl_dirs_file)
yml_paths = wic.main.get_yml_paths(get_args().yml_dirs_file)
yaml_stems = wic.utils.flatten([list(p) for p in yml_paths.values()])
schema_store: Dict[str, Json] = {}
validator = wic.schemas.wic_schema.get_validator(tools_cwl, yaml_stems, schema_store, write_to_disk=True)

yml_paths_tuples = [(yml_path_str, yml_path)
    for yml_namespace, yml_paths_dict in yml_paths.items()
    for yml_path_str, yml_path in yml_paths_dict.items()]

for yml_path_str, yml_path in yml_paths_tuples:
    schema = wic.schemas.wic_schema.compile_workflow_generate_schema(yml_path_str, yml_path,
                                                            tools_cwl, yml_paths, validator)
    # overwrite placeholders in schema_store. See comment in get_validator()
    schema_store[schema['$id']] = schema

validator = wic.schemas.wic_schema.get_validator(tools_cwl, yaml_stems, schema_store, write_to_disk=True)


def wic_yaml_filter_blank_steps(yml: Yaml) -> bool:
    """Filters out Yaml instances with blank steps.

    Args:
        yml (Yaml): A randomly generated Yaml instance.

    Returns:
        bool: True if there are no blank steps.
    """
    steps = yml.get('steps', [])
    return not (steps == [] or any([step == {} for step in steps]))


def wic_yaml_filter_backends_or_steps(yml: Yaml) -> bool:
    """Filters out Yaml instances with no backends and no steps.

    Args:
        yml (Yaml): A randomly generated Yaml instance.

    Returns:
        bool: True if there is either backends or steps (or both).
    """
    return ('backends' in yml or 'steps' in yml)

time_initial = time.time()

wic_schema = wic.schemas.wic_schema.wic_main_schema(tools_cwl, yaml_stems, schema_store, hypothesis=True)
# NOTE: The CLI version of mypy (0.991) and the VSCode version of mypy
# disagree on the following line. According to CLI mypy:
# "Argument 1 to "filter" of "SearchStrategy" has incompatible type"
# If you add a "type: ignore" comment CLI mypy passes, but then VSCode mypy says:
# "Unused "type: ignore" comment"
# Apparently github CI mypy (curiously, also 0.991) agrees with VSCode mypy, so
# do not add a "type: ignore" comment!
wic_strategy = hj.from_schema(wic_schema).filter(wic_yaml_filter_blank_steps).filter(wic_yaml_filter_backends_or_steps)

time_final = time.time()
print(f'from_schema time: {round(time_final - time_initial, 4)} seconds')
print()
