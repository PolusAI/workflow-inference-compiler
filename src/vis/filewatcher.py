import glob
import os
from pathlib import Path
from typing import Dict, List


# NOTE: You should be very careful when using file watchers! Most libraries
# (watchdog, watchfiles, etc) will use operating system / platform-specific
# APIs to check for changes (for performance reasons). However, this can cause
# problems in some cases, specifically for network drives.
# See https://stackoverflow.com/questions/45441623/using-watchdog-of-python-to-monitoring-afp-shared-folder-from-linux
# I'm 99% sure the same problem happens with Docker containers. Either way, the
# solution is to use polling. However, for unknown reasons, simply replacing
# Observer with PollingObserver doesn't seem to work! So we are forced to write
# our own basic file watcher using glob.


def file_watcher_glob(cachedir_path: Path, patterns: List[str], prev_files: Dict[str, float]) -> Dict[str, float]:
    """Determines whether files (specified by the given glob patterns) have been either recently created or modified.\n
    Note that this is a workaround due to an issue with using standard file-watching libraries.

    Args:
        cachedir_path (Path): The --cachedir directory of the main workflow.
        patterns (List[str]): The glob patterns which specifies the files to be watched.
        prev_files (Dict[str, float]): This should be the return value from the previous function call.

    Returns:
        Dict[str, float]: A dictionary containing the filepaths and last modification times.
    """
    changed_files = {}
    file_patterns = [str(cachedir_path / f'**/{pattern}') for pattern in patterns]
    file_paths_all = [glob.glob(fp, recursive=True) for fp in file_patterns]
    file_paths_flat = [path for fps in file_paths_all for path in fps]
    for file in file_paths_flat:
        mtime = os.path.getmtime(file)
        if file not in prev_files:
            # created
            changed_files[file] = mtime
        elif mtime > prev_files[file]:
            # modified
            changed_files[file] = mtime
    return changed_files
