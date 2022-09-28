from pathlib import Path
import time
from typing import Dict, List

from IPython.display import display
import mdtraj
import nglview as nv

from . import filewatcher


def main(num_iterations: int, cachedir_path: str = 'cachedir', file_patterns: List[str] = ['*.trr', '*.pdb']) -> None:
    """This function watches cachedir_path for file_patterns and updates an NGLWidget, upto num_iterations times.

    Args:
        num_iterations (int): The number of iterations, to guarantee termination.
        cachedir_path (str, optional): The cwltool --cachedir directory. Defaults to 'cachedir'.
        file_patterns (List[str], optional): The coordinate and topology file patterns. Defaults to ['*.trr', '*.pdb'].
    """

    prev_files: Dict[str, float] = {}
    coords_file = ''
    top_file = ''
    nglwidget = nv.NGLWidget()
    nglwidget._set_size('100%', '800px') # pylint: disable=protected-access
    display(nglwidget)
    for t in range(num_iterations):
        try:
            time.sleep(1.0)  # Wait at least 1 second so we don't just spin.
            # Use our own polling file watcher, see comments in filewatcher.py
            changed_files = filewatcher.file_watcher_glob(Path(cachedir_path), file_patterns, prev_files)
            changed_files_list = list(changed_files.items())
            changed_files_list.sort(key=lambda x: x[1])
            if len(changed_files_list) == 0:
                continue
            for file, time_ in changed_files_list:
                #print(file)
                if Path(file).suffix == '.trr':
                    coords_file = file
                if Path(file).suffix == '.pdb':
                    top_file = file
            prev_files = {**prev_files, **changed_files}

            #print('coords_file', coords_file)
            #print('top_file', top_file)
            if top_file != '':
                traj: mdtraj.Trajectory
                if coords_file != '':
                    traj = mdtraj.load(coords_file, top=top_file)
                else:
                    traj = mdtraj.load(top_file) # pdb files implicitly contain topology info
                #print(traj)
                #num_frames = 10
                #traj = traj[-num_frames:] # Slice the most recent num_frames
                traj = traj[-1] # Use -1 for most recent frame only

                if nglwidget.n_components == 0: # First time
                    component = nglwidget.add_trajectory(traj)
                else:
                    if nglwidget.max_frame + 1 != traj.n_frames: # max_frame starts from 0
                        # Removing and adding a new trajectory allows increasing
                        # the number of frames. However, this causes the UI to 'blink'
                        nglwidget.remove_component(component)
                        component = nglwidget.add_trajectory(traj)
                    else:
                        # If the number of frames stays the same
                        # (i.e. we're just using the last frame or last n frames)
                        # then we can update the coordinates and avoid the blinking.
                        trajtraj = nv.adaptor.MDTrajTrajectory(traj)
                        for i in range(traj.n_frames):
                            nglwidget.set_coordinates({i: trajtraj.get_coordinates(i)})
        except (IndexError, RuntimeError, OSError, AssertionError, ValueError) as e:
            # TODO: Figure out what is causing these exceptions and determine if anything needs fixed.
            # Some/most of these errors are likely due to race conditions writing/reading to/from disk.
            runtimeerror = 'TRR read error: Float'
            indexerror = 'index 0 is out of bounds for axis 0 with size 0'
            oserror1 = 'Malformed TRR file. Number of atoms <= 0. Are you sure this is a valid GROMACS TRR file?'
            oserror2 = 'No such file:'
            assertionerror = 'assert len(name) == 4'
            ve1 = 'The topology and the trajectory files might not contain the same atoms'
            ve2 = 'xyz must be shape'
            ve3 = 'could not convert string to float:'
            if isinstance(e, RuntimeError) and runtimeerror in str(e):
                pass
            elif isinstance(e, IndexError) and indexerror in str(e):
                pass
            elif isinstance(e, OSError) and (oserror1 in str(e) or oserror2 in str(e)):
                pass
            elif isinstance(e, AssertionError): # and assertionerror in str(e):
                pass
            elif isinstance(e, ValueError) and (ve1 in str(e) or ve2 in str(e) or ve3 in str(e)):
                pass
            else:
                print(str(e))
                raise
    print('done')
