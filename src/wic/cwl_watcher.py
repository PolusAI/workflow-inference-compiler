import glob
import os
from pathlib import Path
import subprocess as sub
import sys
import time
from typing import Any, Dict, List
from unittest.mock import patch

import graphviz
import networkx as nx
import json

"""from watchdog.observers import Observer
from watchdog.observers.polling import PollingObserver
from watchdog.events import FileSystemEvent, PatternMatchingEventHandler"""

from .wic_types import GraphReps, Tools, YamlTree, Json
from .main import get_tools_cwl

from . import cli, compiler, utils


def absolute_paths(json: Json, path: str) -> Json:
    new_json: Json = {}
    for key, val in json.items():
        new_val = val
        # TODO: Generalize beyond the gmx_energy use case
        if key == 'input_energy_path':
            new_val = path + '/' + val
        new_json[key] = new_val
    return new_json


def rerun_cwltool(directory: Path, cachedir_path: str, cwl_tool: str, args_vals: Json, tools_cwl: Tools) -> None:
    try:
        # Make paths in arguments absolute w.r.t the given directory. See below.
        args_vals_new = absolute_paths(args_vals, str(directory))

        # Construct a single-step workflow and add its arguments
        yml = {'steps': [{cwl_tool: {'in': args_vals_new}}]}
        #print('yml')
        #print(yml)
        #import yaml
        #print(yaml.dump(yml))

        # Measure compile time
        time_initial = time.time()

        # Setup dummy args
        testargs = ['wic', '--yaml', '', '--cwl_output_intermediate_files', 'True']  # ignore --yaml
        # For now, we need to enable --cwl_output_intermediate_files. See comment in compiler.py
        with patch.object(sys, 'argv', testargs):
            args = cli.parser.parse_args()

        yaml_path = f'{cwl_tool}_only.yml'
        yaml_tree = YamlTree(yaml_path, yml)
        subgraph = GraphReps(graphviz.Digraph(name=yaml_path), nx.DiGraph())
        compiler_info = compiler.compile_workflow(yaml_tree, args, [], [subgraph], {}, {}, tools_cwl, True, relative_run_path=False)
        rose_tree = compiler_info.rose
        working_dir = Path('.') # Use a new working directory.
        # Can also use `directory` at the risk of overwriting other files.
        utils.write_to_disk(rose_tree, working_dir, relative_run_path=False)

        time_final = time.time()
        print(f'compile time for {cwl_tool}: {round(time_final - time_initial, 4)} seconds')

        # NOTE: Since we are running cwltool 'within' cwltool, the inner
        # cwltool command will get run from working_dir, but then cwl_tool
        # will run within some other hashed directory in .../cachedir/
        # The solution is to modify the input paths above to be absolute.
        # The easiest way to do this for now is recompiling. This adds a few
        # seconds, but most of the time will be CWL validation and runtime.
        # Alternatively, we could try to compile once in main() and then
        # make the paths absolute in f'{cwl_tool}_only_inputs.yml' here.
        cmd: List[str] = ['cwltool', '--cachedir', cachedir_path, f'{cwl_tool}_only.cwl', f'{cwl_tool}_only_inputs.yml']
        #proc = sub.run(self.cmd, cwd=working_dir)
        #cmd = self.cmd
        print('Running', cmd)
        proc = sub.run(cmd, cwd=working_dir)
        print('inner cwltool completed')
        # Don't check the return code because the file may not exist yet, or
        # because speculative execution may fail for any number of reasons.
        # proc.check_returncode()
    except FileNotFoundError as e:
        # The file may not exist yet.
        print(e)


# NOTE: You should be very careful when using file watchers! Most libraries
# (watchdog, watchfiles, etc) will use operating system / platform-specific
# APIs to check for changes (for performance reasons). However, this can cause
# problems in some cases, specifically for network drives.
# See https://stackoverflow.com/questions/45441623/using-watchdog-of-python-to-monitoring-afp-shared-folder-from-linux
# I'm 99% sure the same problem happens with Docker containers. Either way, the
# solution is to use polling. However, for unknown reasons, simply replacing
# Observer with PollingObserver doesn't seem to work! So we are forced to write
# our own basic file watcher using glob.


"""class SubprocessHandler(PatternMatchingEventHandler):

    def __init__(self, cmd: List[str], cachedir_path: str, cwl_tool: str, args_vals: Json, tools_cwl: Tools, **kwargs: Any) -> None:
        self.cmd = cmd
        self.lock = False
        self.cachedir_path = cachedir_path
        self.cwl_tool = cwl_tool
        self.args_vals = args_vals
        self.tools_cwl = tools_cwl
        super().__init__(**kwargs)

    def on_any_event(self, event: FileSystemEvent) -> None:
        # Use a lock to prevent us from DOS'ing ourselves
        global lock
        if event.event_type == 'modified' and not lock:
            directory = Path(event._src_path).parent
            print('directory', directory)
            #self.lock = True
            print(event)
            rerun_cwltool(directory, self.cachedir_path, self.cwl_tool, self.args_vals, self.tools_cwl)
            #self.lock = False"""


def file_watcher_glob(dir: Path, pattern: str, prev_files: Dict[str, float]) -> Dict[str, float]:
    changed_files = {}
    file_pattern = str(dir / f'**/{pattern}')
    file_paths = glob.glob(file_pattern, recursive=True)
    for file in file_paths:
        time = os.path.getmtime(file)
        if file not in prev_files:
            # created
            changed_files[file] = time
        elif time > prev_files[file]:
            # modified
            changed_files[file] = time
    return changed_files


def main() -> None:
    print('cwl_watcher sys.argv', sys.argv)
    # TODO: check that 1,3,5, are --cachedir_path, --file_pattern, --cwl_tool
    # or switch to argparse
    cachedir_path = sys.argv[2]
    file_pattern = sys.argv[4].strip()
    cwl_tool = sys.argv[6]
    max_times = int(sys.argv[8])

    # Create an empty 'logfile' so that cwl_watcher.cwl succeeds.
    # TODO: Maybe capture cwl_tool stdout/stderr and redirect to this logfile.
    logfile = Path(f'{cwl_tool}_only.log')
    logfile.touch()

    # Parse config into CWL input args
    config = sys.argv[10]
    args_vals = json.loads(config)

    # This really needs to be args.cwl_dir, where args comes from the original
    # command line, i.e. we can't just use dummy args.
    cwl_dir = Path(cachedir_path).parent
    tools_cwl = get_tools_cwl(cwl_dir)

    """# Make paths in arguments absolute w.r.t the given directory. See below.
    args_vals_new = absolute_paths(args_vals, str(directory))

    # Construct a single-step workflow and add its arguments
    yml = {'steps': [{cwl_tool: {'in': args_vals_new}}]}
    #print('yml')
    #print(yml)
    #print(yaml.dump(yml))

    # Setup dummy args
    testargs = ['wic', '--yaml', '', '--cwl_output_intermediate_files', 'True']  # ignore --yaml
    # For now, we need to enable --cwl_output_intermediate_files. See comment in compiler.py
    with patch.object(sys, 'argv', testargs):
        args = cli.parser.parse_args()

    yaml_path = f'{cwl_tool}_only.yml'
    yaml_tree = YamlTree(yaml_path, yml)
    subgraph = GraphReps(graphviz.Digraph(name=yaml_path), nx.DiGraph())
    compiler_info = compiler.compile_workflow(yaml_tree, args, [], [subgraph], {}, {}, tools_cwl, True, relative_run_path=False)
    rose_tree = compiler_info.rose
    utils.write_to_disk(rose_tree, Path('.'), relative_run_path=False)"""

    cachedir_hash_path = Path('.').absolute()
    print('cachedir_hash_path', cachedir_hash_path)

    """cmd: List[str] = ['cwltool', '--cachedir', cachedir_path, f'{cachedir_hash_path}/{cwl_tool}_only.cwl', f'{cachedir_hash_path}/{cwl_tool}_only_inputs.yml']
    event_handler = SubprocessHandler(cmd, cachedir_path, cwl_tool, args_vals, tools_cwl, patterns=[file_pattern])
    observer = PollingObserver()  # This does not work!
    observer.schedule(event_handler, cachedir_path, recursive=True)
    observer.start()"""

    # Specify a maximum number of iterations to guarantee termination.
    # Total runtime will be (sleep time + compile time + run time) * max_iters
    # For now, there is no way to estimate max_iters such that polling will end
    # around the same time as the original workflow step.
    # TODO: Generate a file when the original workflow step finishes, and look
    # for that file here to terminate. Keep max_iters just in case.
    iter = 0
    prev_files: Dict[str, float] = {}
    try:
        while iter < max_times:
            # Use our own polling file watcher, see above.
            changed_files = file_watcher_glob(cachedir_hash_path.parent, '*.edr', prev_files)
            for file in changed_files:
                if file_pattern[1:] in file:
                    print(file)
                    rerun_cwltool(Path(file).parent, cachedir_path, cwl_tool , args_vals, tools_cwl)
            prev_files = {**prev_files, **changed_files}

            time.sleep(1.0) # Wait at least 1 second so we don't just spin.
            iter += 1
    except KeyboardInterrupt:
        pass
    #observer.stop()
    #observer.join()

    failed = False # Your analysis goes here
    if failed:
        print(f'{cwl_tool} failed!')
        sys.exit(1)


if __name__ == "__main__":
    main()