import argparse
import copy
import json
from pathlib import Path
from typing import Any, List, Tuple

import yaml

from . import auto_gen_header
from .wic_types import (Namespaces, NodeData, RoseTree, Yaml, ExplicitEdgeCalls, Json)


def read_lines_pairs(filename: Path) -> List[Tuple[str, str]]:
    """Reads a whitespace-delimited file containing two paired entries per line (i.e. a serialized Dict).

    Args:
        filename (Path): The full path of the file to be read.

    Raises:
        Exception: If any non-blank, non-comment lines do not contain exactly two entries.

    Returns:
        List[Tuple[str, str]]: The file contents, with blank lines and comments removed.
    """
    with open(filename, mode='r', encoding='utf-8') as f:
        lines = []
        for line in f.readlines():
            if line.strip() == '':  # Skip blank lines
                continue
            if line.startswith('#'):  # Skip comment lines
                continue
            l_s = line.split()
            if not len(l_s) == 2:
                print(line)
                raise Exception("Error! Line must contain exactly two entries!")
            lines.append((l_s[0], l_s[1]))
    return lines


# snakeyaml (a cromwell dependency) refuses to parse yaml files with more than
# 50 anchors/aliases to prevent Billion Laughs attacks.
# See https://en.wikipedia.org/wiki/Billion_laughs_attack
# Solution: Inline the contents of the aliases into the anchors.
# See https://ttl255.com/yaml-anchors-and-aliases-and-how-to-disable-them/#override


class NoAliasDumper(yaml.SafeDumper):
    def ignore_aliases(self, data: Any) -> bool:
        return True


def write_to_disk(rose_tree: RoseTree, path: Path, relative_run_path: bool) -> None:
    """Writes the compiled CWL files and their associated yml inputs files to disk.

    NOTE: Only the yml input file associated with the root workflow is
    guaranteed to have all inputs. In other words, subworkflows will all have
    valid CWL files, but may not be executable due to 'missing' inputs.

    Args:
        rose_tree (RoseTree): The data associated with compiled subworkflows
        path (Path): The directory in which to write the files
        relative_run_path (bool): Controls whether to use subdirectories or just one directory.
    """
    node_data: NodeData = rose_tree.data
    namespaces = node_data.namespaces
    yaml_stem = node_data.name
    cwl_tree = node_data.compiled_cwl
    yaml_inputs = node_data.workflow_inputs_file

    # NOTE: As part of the scatter feature we introduced the use of 'source',
    # but in some cases (biobb 'config' tag) it is not being removed correctly
    # in the compiler, so as a last resort remove it here.
    yaml_inputs_no_source = {}
    for key, val in yaml_inputs.items():
        try:
            if isinstance(val, str):
                val_dict = json.loads(val)
                if 'source' in val_dict:
                    val = val_dict['source']
        except Exception as e:
            pass
        yaml_inputs_no_source[key] = val

    path.mkdir(parents=True, exist_ok=True)
    if relative_run_path:
        filename_cwl = f'{yaml_stem}.cwl'
        filename_yml = f'{yaml_stem}_inputs.yml'
    else:
        filename_cwl = '___'.join(namespaces + [f'{yaml_stem}.cwl'])
        filename_yml = '___'.join(namespaces + [f'{yaml_stem}_inputs.yml'])

    # Dump the compiled CWL file contents to disk.
    # Use sort_keys=False to preserve the order of the steps.
    yaml_content = yaml.dump(cwl_tree, sort_keys=False, line_break='\n', indent=2, Dumper=NoAliasDumper)
    with open(path / filename_cwl, mode='w', encoding='utf-8') as w:
        w.write('#!/usr/bin/env cwl-runner\n')
        w.write(auto_gen_header)
        w.write(''.join(yaml_content))

    yaml_content = yaml.dump(yaml_inputs_no_source, sort_keys=False, line_break='\n', indent=2, Dumper=NoAliasDumper)
    with open(path / filename_yml, mode='w', encoding='utf-8') as inp:
        inp.write(auto_gen_header)
        inp.write(yaml_content)

    for sub_rose_tree in rose_tree.sub_trees:
        subpath = path
        if relative_run_path:
            sub_node_data: NodeData = sub_rose_tree.data
            sub_step_name = sub_node_data.namespaces[-1]
            subpath = path / sub_step_name
        write_to_disk(sub_rose_tree, subpath, relative_run_path)


def write_config_to_disk(config: Json, config_file: Path) -> None:
    """Writes config json object to config_file

    Args:
        config (Json): The json object that is to be written to disk
        config_file (Path): The file path where it is to be written
    """
    config_dir = Path(config_file).parent
    # make the full path if it doesn't exist
    config_dir.mkdir(parents=True, exist_ok=True)
    with open(config_file, 'w', encoding='utf-8') as f:
        json.dump(config, f)


def read_config_from_disk(config_file: Path) -> Json:
    """Returns the config json object from config_file with absolute paths

    Args:
        config_file (Path): The path of json file where it is to be read from

    Returns:
        Json: The config json object with absolute filepaths
    """
    config: Json = {}
    # config_file can contain absolute or relative paths
    with open(config_file, 'r', encoding='utf-8') as f:
        config = json.load(f)
    conf_tags = ['search_paths_cwl', 'search_paths_yml']
    for tag in conf_tags:
        config[tag] = get_absolute_paths(config[tag])
    return config


def get_default_config() -> Json:
    """Returns the default config with absolute paths

    Returns:
        Json: The config json object with absolute filepaths
    """
    src_dir = Path(__file__).parent
    conf_tags = ['search_paths_cwl', 'search_paths_yml']
    default_config: Json = {}
    # config.json can contain absolute or relative paths
    default_config = read_config_from_disk(src_dir/'config.json')
    for tag in conf_tags:
        default_config[tag] = get_absolute_paths(default_config[tag])
    return default_config


def get_absolute_paths(sub_config: Json) -> Json:
    """Makes the paths within the dirs_file file absolute and write them into sub_config object.

    Args:
        sub_config (dict): The json (sub)object where filepaths are stored

    Returns:
        Json: The json (sub)object with absolute filepaths
    """
    abs_sub_config = copy.deepcopy(sub_config)
    for ns in abs_sub_config:
        abs_paths = [str(Path(path).absolute()) for path in abs_sub_config[ns]]
        abs_sub_config[ns] = abs_paths
    return abs_sub_config


def write_absolute_yaml_tags(args: argparse.Namespace, in_dict_in: Yaml, namespaces: Namespaces,
                             step_name_i: str, explicit_edge_calls_copy: ExplicitEdgeCalls) -> None:
    """cwl_watcher requires all paths to be absolute.

    Args:
        args (argparse.Namespace): The command line arguments
        in_dict_in (Yaml): The in: subtag of a cwl_watcher: tag. (Mutates in_dict_in)
        namespaces (Namespaces): Specifies the path in the yml AST to the current subworkflow
        step_name_i (str): The name of the current workflow step
        explicit_edge_calls_copy (ExplicitEdgeCalls): Stores the (path, value) of the explicit edge call sites
    """

    # cachedir_path needs to be an absolute path, but for reproducibility
    # we don't want users' home directories in the yml files.
    cachedir_path = Path(args.cachedir).absolute()
    # print('setting cachedir_path to', cachedir_path)
    in_dict_in['root_workflow_yml_path'] = str(Path(args.yaml).parent.absolute())

    in_dict_in['cachedir_path'] = str(cachedir_path)
    in_dict_in['homedir'] = args.homedir

    # Add a 'dummy' values to explicit_edge_calls, because
    # that determines sub_args_provided when the recursion returns.
    arg_keys_ = ['root_workflow_yml_path', 'cachedir_path', 'homedir']
    for arg_key_ in arg_keys_:
        in_name_ = f'{step_name_i}___{arg_key_}'  # {step_name_i}_input___{arg_key}
        explicit_edge_calls_copy.update({in_name_: (namespaces + [step_name_i], arg_key_)})
