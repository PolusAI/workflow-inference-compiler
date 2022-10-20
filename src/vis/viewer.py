import json
from pathlib import Path
import time
import threading
from typing import Dict, List

from IPython.display import display
import mdtraj
import nglview as nv

from ipytree import Tree, Node
from ipywidgets import HBox

from wic import utils
from . import filewatcher


def make_ipytree_nodes(prov_tree: Dict, rootdir: str) -> List[Node]:
    """Recursively applies the ipytree Node constructor to each element of the provenance output files tree.

    Args:
        prov_tree (Dict): This should be the output of utils.provenance_list_to_tree(files)
        rootdir (str): The directory in which to search for outdir/

    Returns:
        List[Node]: A list of Nodes (where each Node contains children Nodes, etc)
    """
    nodes = []
    for key, val in prov_tree.items():
        if isinstance(val, Dict):
            children = make_ipytree_nodes(val, rootdir)
        if isinstance(val, List):
            children = []
            dests = set()
            for (location, parentdirs, basename) in val:
                dest = rootdir + 'outdir/' + parentdirs +  '/' + basename
                idx = 2
                while dest in dests:
                    stem = Path(basename).stem
                    suffix = Path(basename).suffix
                    basename_ = stem + f'_{idx}' + suffix
                    dest = rootdir + 'outdir/' + parentdirs +  '/' + basename_
                    idx += 1
                if idx > 2: # If the while loop ran at least once
                    basename = basename_
                dests.add(dest)
                child = Node(basename)
                # NOTE: We need to store a unique id in each node so we can
                # take the appropriate action when the user clicks. Storing
                # the id in the name attribute would make the UI look terrible.
                # However, this is python, so we can just pretend an id
                # attribute exists and use it anyway. Is this unsafe? Probably!
                child.id = parentdirs
                children.append(child)
        node = Node(key, children)
        node.id = '' # See above comment.
        node.opened = False
        nodes.append(node)
    return nodes


def tree_viewer(rootdir: str = '../../') -> HBox:
    """Creates a file browser that will display molecular files.

    Args:
        rootdir (str, optional): The directory in which to search for\n
        provenance/workflow/primary-output.json'. Defaults to '../../'.

    Returns:
        HBox: An ipywidget which contains ipytree on the left and nglview on the right.
    """
    nglwidget = nv.NGLWidget()
    nglwidget._set_size('100%', '800px') # pylint: disable=protected-access

    output_json_file = Path(rootdir + 'provenance/workflow/primary-output.json')
    if not output_json_file.exists():
        print(f'Error! {output_json_file.absolute()} does not exist!')
        print('Did your workflow finish executing?')
        return HBox()

    files = utils.parse_provenance_output_files(output_json_file)
    with open(output_json_file, mode='r', encoding='utf-8') as f:
        output_dict = json.loads(f.read())
    prov_tree = utils.provenance_list_to_tree(files)
    #import yaml
    #print(yaml.dump(prov_tree))
    children = make_ipytree_nodes(prov_tree, rootdir)

    rootnode = Node("Workflow", children)
    tree  = Tree(nodes=[rootnode], multiple_selection=False)
    components: List = []

    def clear_components() -> None:
        for component in components:
            nglwidget.remove_component(component)
        components.clear()

    def on_selected_change(change: Dict) -> None:
        #print('change[new]', change['new'])
        parentdirs = change['new'][0].id
        basename = change['new'][0].name
        #print('parentdirs:', parentdirs)
        if parentdirs != '':
            filepath = rootdir + 'outdir/' + parentdirs + '/' + basename
            #print(filepath)
            if Path(filepath).exists():
                mdtraj_exts = [".pdb", ".pdb.gz", ".h5", ".lh5", ".prmtop", ".parm7", ".prm7",
                               ".psf", ".mol2", ".hoomdxml", ".gro", ".arc", ".hdf5", ".gsd"]
                # NOTE: mdtraj does not support .sdf and .pdbqt
                ngl_exts = [".mmcif", ".cif", ".mcif", ".pdb", ".pdbqt", ".ent",
                            ".pqr", ".gro", ".sdf", ".sd", ".mol2", ".mmtf"]
                if Path(filepath).suffix in mdtraj_exts:
                    clear_components()
                    traj = mdtraj.load(filepath)
                    component = nglwidget.add_trajectory(traj)
                    components.append(component)
                elif Path(filepath).suffix in ngl_exts:
                    clear_components()
                    component = nglwidget.add_component(filepath)
                    components.append(component)

    tree.observe(on_selected_change, names='selected_nodes')
    return HBox([tree, nglwidget])


def realtime_viewer(num_iterations: int, cachedir_path: str = '../../cachedir', file_patterns: List[str] = ['*.trr', '*.pdb']) -> None:
    """This function watches cachedir_path for file_patterns and updates an NGLWidget, upto num_iterations times.

    Args:
        num_iterations (int): The number of iterations, to guarantee termination.
        cachedir_path (str, optional): The cwltool --cachedir directory. Defaults to 'cachedir'.
        file_patterns (List[str], optional): The coordinate and topology file patterns. Defaults to ['*.trr', '*.pdb'].
    """
    # Just like matplotlib, you can't run a calculation in the GUI event loop thread
    # or else the GUI will not redraw. However, once the realtime_viewer() thread finishes
    # there is no easy way to interrupt another thread (i.e. no Ctrl-C), so simply
    # use a fixed number of iterations.
    thread = threading.Thread(target=realtime_viewer_body, args=(num_iterations, cachedir_path, file_patterns))
    thread.start()


def realtime_viewer_body(num_iterations: int, cachedir_path: str, file_patterns: List[str]) -> None:
    """This function watches cachedir_path for file_patterns and updates an NGLWidget, upto num_iterations times.

    Args:
        num_iterations (int): The number of iterations, to guarantee termination.
        cachedir_path (str, optional): The cwltool --cachedir directory.
        file_patterns (List[str], optional): The coordinate and topology file patterns.
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
