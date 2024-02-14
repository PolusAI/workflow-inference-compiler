import argparse
import sys
import time
from typing import Dict, List
from unittest.mock import patch
from pathlib import Path

from hypothesis.strategies import SearchStrategy
import hypothesis_jsonschema as hj

import wic
import wic.cli
import wic.plugins
import wic.schemas
import wic.schemas.wic_schema
import wic.utils
from wic.wic_types import Json, Yaml


def get_args(yaml_path: str = '', cwl_inline_subworkflows: bool = False) -> argparse.Namespace:
    """This is used to get mock command line arguments.

    Returns:
        argparse.Namespace: The mocked command line arguments
    """
    testargs = ['wic', '--yaml', yaml_path, '--cwl_output_intermediate_files', 'True']  # ignore --yaml
    if cwl_inline_subworkflows:
        testargs.append("--cwl_inline_subworkflows")
    # For now, we need to enable --cwl_output_intermediate_files. See comment in compiler.py
    with patch.object(sys, 'argv', testargs):
        args: argparse.Namespace = wic.cli.parser.parse_args()
    return args


def get_appended_args(yaml_path: str = '', suppliedargs: List[str] = []) -> argparse.Namespace:
    """This is used to get mock command line arguments, default + suppled args

    Returns:
        argparse.Namespace: The mocked command line arguments
    """
    defaultargs = ['wic', '--yaml', yaml_path, '--cwl_output_intermediate_files', 'True']  # ignore --yaml
    testargs = defaultargs + suppliedargs
    # For now, we need to enable --cwl_output_intermediate_files. See comment in compiler.py
    with patch.object(sys, 'argv', testargs):
        args: argparse.Namespace = wic.cli.parser.parse_args()
    return args


args = get_args()
tools_cwl = wic.plugins.get_tools_cwl(args.homedir, quiet=args.quiet)
yml_paths = wic.plugins.get_yml_paths(args.homedir)
yaml_stems = wic.utils.flatten([list(p) for p in yml_paths.values()])
schema_store: Dict[str, Json] = {}
validator = wic.schemas.wic_schema.get_validator(tools_cwl, yaml_stems, schema_store, write_to_disk=True)

yml_paths_tuples = [(yml_path_str, yml_path)
                    for yml_namespace, yml_paths_dict in yml_paths.items()
                    for yml_path_str, yml_path in yml_paths_dict.items()]

yml_paths_partial_failure = [('examples/test_rand_fail.yml', Path('examples/test_rand_fail.yml'))]

for yml_path_str, yml_path in yml_paths_tuples:
    schema = wic.schemas.wic_schema.compile_workflow_generate_schema(args.homedir, yml_path_str, yml_path,
                                                                     tools_cwl, yml_paths, validator,
                                                                     args.ignore_validation_errors)
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
wic_strategy: SearchStrategy = hj.from_schema(wic_schema)
# NOTE: The CLI version of mypy and the VSCode version of mypy disagree on the
# following line. The "type: ignore" comment is NOT unused.
wic_strategy = wic_strategy.filter(wic_yaml_filter_blank_steps)
wic_strategy = wic_strategy.filter(wic_yaml_filter_backends_or_steps)

time_final = time.time()
print(f'from_schema time: {round(time_final - time_initial, 4)} seconds')
print()
