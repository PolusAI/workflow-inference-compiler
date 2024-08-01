import json
from pathlib import Path
from shutil import copytree, ignore_patterns
from setuptools import setup, find_packages
from setuptools.command.install import install
import versioneer
# from sophios import _version, cli, input_output as io


class PostInstallCommand(install):
    """Post-installation for installation mode."""

    def run(self) -> None:
        install.run(self)
        # Custom post-installation actions go here
        # Never overwrite user config
        if not (Path.home() / 'wic').exists():
            # config_path
            config_path = Path.home() / 'wic'
            cwl_path = config_path / 'cwl_adapters'
            wic_path = config_path / 'examples'
            global_config_file = config_path / 'global_config.json'
            extlist = ['*.png', '*.md', '*.rst', '*.pyc', '__pycache__', '*.json']
            adapters_dir = Path(__file__).parent / 'cwl_adapters'
            examples_dir = Path(__file__).parent / 'docs' / 'tutorials'
            config_file = Path(__file__).parent / 'src' / 'sophios' / 'config_basic.json'
            config = {}
            # config_file can contain absolute or relative paths
            with open(config_file, 'r', encoding='utf-8') as f:
                config = json.load(f)
            conf_tags = ['search_paths_cwl', 'search_paths_wic']
            for tag in conf_tags:
                for ns in config[tag]:
                    abs_paths = [str(Path.home() / path) for path in config[tag][ns]]
                    config[tag][ns] = abs_paths
            # copy to config_path
            copytree(adapters_dir, cwl_path)
            copytree(examples_dir, wic_path, ignore=ignore_patterns(*extlist))
            with open(global_config_file, 'w', encoding='utf-8') as f:
                json.dump(config, f)


setup(
    name='sophios',
    version=versioneer.get_version(),
    packages=find_packages(),
    cmdclass={
        'install': PostInstallCommand,
    },
)
