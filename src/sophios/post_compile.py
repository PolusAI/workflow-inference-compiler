import argparse
from pathlib import Path
import subprocess as sub
from typing import Dict, Union
from . import plugins
from .wic_types import RoseTree


def find_output_dirs(data: Union[RoseTree, Dict, list]) -> list:
    """
    Recursively searches through a nested structure and finds all dictionaries 
    that contain the key 'location', and a key 'class' with a value of 'Directory'.

    Args:
        data (any): The data to search through, which can be a dictionary, list, 
        or any other structure.

    Returns:
        list: A list of location values.
    """
    results = []
    if isinstance(data, Dict):
        if "class" in data and data["class"] == "Directory" and "location" in data:
            results.append(data["location"])
        for value in data.values():
            results.extend(find_output_dirs(value))
    elif isinstance(data, list):
        for item in data:
            results.extend(find_output_dirs(item))

    return results


def create_output_dirs(output_dicts: list, basepath: str = 'autogenerated') -> None:
    """
    Creates all the directories that are needed for the outputs of a workflow.
    """
    for output_dict in output_dicts:
        ldir = Path(output_dict)
        if not ldir.is_absolute():
            ldir = Path(basepath) / ldir
        ldir.mkdir(parents=True, exist_ok=True)


def find_and_create_output_dirs(rose_tree: RoseTree, basepath: str = 'autogenerated') -> None:
    """
    Finds all output directories in the workflow and creates them.
    """
    output_dirs = find_output_dirs(rose_tree.data.workflow_inputs_file)
    create_output_dirs(output_dirs, basepath)


def cwl_docker_extract(args: argparse.Namespace, file_name: str) -> None:
    """Helper function to do the cwl_docker_extract"""
    # cwl-docker-extract recursively `docker pull`s all images in all subworkflows.
    # This is important because cwltool only uses `docker run` when executing
    # workflows, and if there is a local image available,
    # `docker run` will NOT query the remote repository for the latest image!
    # cwltool has a --force-docker-pull option, but this may cause multiple pulls in parallel.
    if args.container_engine == 'singularity':
        cmd = ['cwl-docker-extract', '-s', '--dir',
               f'{args.singularity_pull_dir}', f'autogenerated/{file_name}.cwl']
    else:
        cmd = ['cwl-docker-extract', '--force-download', f'autogenerated/{file_name}.cwl']
    sub.run(cmd, check=True)


def cwl_inline_runtag(args: argparse.Namespace, rose_tree: RoseTree) -> RoseTree:
    """Transform with cwl inline runtag"""
    # this has to happen after at least one write
    # so we can copy from local cwl_dapters in autogenerated/
    if args.cwl_inline_runtag:
        rose_tree = plugins.cwl_update_inline_runtag_rosetree(rose_tree, Path('autogenerated/'), True)
    return rose_tree


def remove_entrypoints(args: argparse.Namespace, rose_tree: RoseTree) -> RoseTree:
    """Remove entry points"""
    if args.docker_remove_entrypoints:
        # Requires root, so guard behind CLI option
        if args.container_engine == 'docker':
            plugins.remove_entrypoints_docker()
        if args.container_engine == 'podman':
            plugins.remove_entrypoints_podman()

        rose_tree = plugins.dockerPull_append_noentrypoint_rosetree(rose_tree)
    return rose_tree
