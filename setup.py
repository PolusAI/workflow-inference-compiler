import json
from pathlib import Path
from shutil import copytree, ignore_patterns
from setuptools import setup, find_packages
from setuptools.command.build_py import build_py as _build_py
import versioneer


class build_py(_build_py):
    """build mode script"""

    def run(self) -> None:
        adapters_dir = Path(__file__).parent / 'cwl_adapters'
        examples_dir = Path(__file__).parent / 'docs' / 'tutorials'
        adapters_target_dir = Path(self.build_lib) / 'sophios/cwl_adapters'
        examples_target_dir = Path(self.build_lib) / 'sophios/examples'
        extlist = ['*.png', '*.md', '*.rst', '*.pyc', '__pycache__', '*.json']
        copytree(adapters_dir, adapters_target_dir, dirs_exist_ok=True)
        copytree(examples_dir, examples_target_dir, ignore=ignore_patterns(*extlist), dirs_exist_ok=True)
        # Never overwrite user config
        if not (Path.home() / 'wic').exists():
            # config_path
            config_path = Path.home() / 'wic'
            # global config file in ~/wic
            global_config_file = config_path / 'global_config.json'

            # setup.py always gets executed in the project_root
            config_file = Path(__file__).parent / 'src' / 'sophios' / 'config.json'
            config = {}
            # config_file can contain absolute or relative paths
            with open(config_file, 'r', encoding='utf-8') as f:
                config = json.load(f)
            conf_tags = ['search_paths_cwl', 'search_paths_wic']
            for tag in conf_tags:
                for ns in config[tag]:
                    abs_paths = [str(Path(path).absolute()) for path in config[tag][ns]]
                    config[tag][ns] = abs_paths
            config_path.mkdir(parents=True, exist_ok=True)
            # write out the config file with paths
            with open(global_config_file, 'w', encoding='utf-8') as f:
                json.dump(config, f)

        # Continue with the standard build process
        super().run()


setup(
    name='sophios',
    version=versioneer.get_version(),
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    include_package_data=True,
    cmdclass={
        'build_py': build_py,
    },
    package_data={
        "sophios": ["cwl_adapters/*.cwl", "examples/*.wic"]
    }
)
