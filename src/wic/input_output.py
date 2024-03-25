import argparse
import logging
import json
from pathlib import Path
import subprocess as sub
from typing import Any, List, Tuple

import yaml

from . import auto_gen_header
from .wic_types import (Namespaces, NodeData, RoseTree, Yaml, ExplicitEdgeCalls)


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


logger_wicad = logging.getLogger("wicautodiscovery")


def copy_config_files(homedir: str) -> None:
    """Copies the following configuration files to ~/wic/\n
    cwl_dirs.txt, yml_dirs.txt, renaming_conventions.txt, inference_rules.txt

    Args:
        homedir (str): The users home directory
    """
    files = ['cwl_dirs.txt', 'yml_dirs.txt', 'renaming_conventions.txt', 'inference_rules.txt']
    src_dir = Path(__file__).parent
    wicdir = Path(homedir) / 'wic'
    wicdir.mkdir(exist_ok=True)

    for file in files:
        if not (wicdir / file).exists():
            logger_wicad.warning(f'Writing {str(wicdir / file)}')
            logger_wicad.warning('Please check this file and make sure that the paths in it are correct.')
            cmd = ['cp', str(src_dir / file), str(wicdir / file)]
            sub.run(cmd, check=True)

    write_absolute_config_files(wicdir / 'cwl_dirs.txt')
    write_absolute_config_files(wicdir / 'yml_dirs.txt')


def write_absolute_config_files(dirs_file: Path) -> None:
    """Makes the paths within the \*_dirs.txt files absolute

    Args:
        dirs_file (Path): The path to the \*_dirs.txt file
        dirs_file_abs (str): The path to the absolute \*_dirs.txt file
    """
    ns_paths = read_lines_pairs(dirs_file)
    pairs_abs = [ns + ' ' + str(Path(path).absolute()) for ns, path in ns_paths]
    with open(dirs_file, mode='w', encoding='utf-8') as f:
        f.write('\n'.join(pairs_abs))


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

    write_absolute_config_files(Path(args.homedir) / 'wic' / 'cwl_dirs.txt')
    write_absolute_config_files(Path(args.homedir) / 'wic' / 'yml_dirs.txt')
