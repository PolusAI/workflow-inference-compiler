import shutil
from pathlib import Path
from setuptools import setup, find_packages
from setuptools.command.build_py import build_py as _build_py
import versioneer

class build_py(_build_py):
    def run(self):
        # Copy cwl adapters to the package directory
        target_dir = Path(self.build_lib) / 'sophios/cwl_adapters'
        adapters_dir = Path(__file__).parent / 'cwl_adapters'
        shutil.copytree(adapters_dir, target_dir, dirs_exist_ok=True)
        
        # Continue with the standard build process
        super().run()

setup(
    name="sophios",
    version=versioneer.get_version(),
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    include_package_data=True,
    cmdclass={'build_py': build_py},
    package_data={
        "sophios": ["cwl_adapters/*.cwl"], 
    },
)
