import glob
from pathlib import Path
import os
import sys
import time
from typing import List, Dict, Tuple
import matplotlib

import numpy as np
from matplotlib import _pylab_helpers # type: ignore
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

import ruptures as rpt

"""from watchdog.observers import Observer
from watchdog.observers.polling import PollingObserver
from watchdog.events import FileSystemEvent, PatternMatchingEventHandler"""


def read_tabular_data(filename: Path) -> List[List[float]]:
    lines = []
    for line in open(filename, 'r').readlines():
        if line.startswith('#'):  # Skip comment lines
            continue
        lines.append([float(x) for x in line.strip().split()])
    return lines


# Declare global variables to pass data between the main thread and the observer thread. (Not ideal)
data_glob: List[Tuple[Path, np.ndarray]] = []
data_glob_changed: bool = False


def store_tabular_data(path: Path, use_stem: bool = True) -> None:
    # Declare global variables locally
    global data_glob
    global data_glob_changed
    floats = read_tabular_data(path)
    
    if floats == []:
        print('Skipping empty file', path)
        return None

    # Check that the array is not ragged; each line must be the same length!
    # I'm not exactly sure why this happens, but it seems like maybe the file
    # contents are not being flushed to disk before getting read back in again.
    # When I manually check the files afterwards, the data is all present.
    lengths = [len(x) for x in floats]
    if not all([length == lengths[0] for length in lengths]):
        print('Warning! Skipping ragged data in', path)
        return None

    data = np.array(floats)
    if use_stem:
        path = Path(path.stem)
    for i in range(len(data_glob)):
        (p, data_old_) = data_glob[i]
        if path == p:
            data_glob[i] = (path, data)
            data_glob_changed = True
            return None
    data_glob.append((path, data))
    data_glob_changed = True
    return None


def initialize_plots(nrows: int, ncols: int) -> Tuple[matplotlib.pyplot.Figure, List[List[matplotlib.pyplot.Axes]]]:
    """Initialize a grid of blank subplots.

    Args:
        nrows (int): the number of rows
        ncols (int): the number of columns

    Returns:
        Tuple[matplotlib.pyplot.Figure, List[List[matplotlib.pyplot.Axes]]]: The main figure and subplots axes
    """
    plt.style.use('dark_background') # type: ignore
    fig, axes2d = plt.subplots(nrows=nrows, ncols=ncols, figsize=(18, 9.5))
    fig.suptitle('')
    fig.tight_layout()
    plt.subplots_adjust(left=0.05, wspace=0.24, bottom=0.05, hspace=0.24)
    plt.show(block=False) # type: ignore
    return (fig, axes2d)


# For now, hardcode labels for the gromacs example.
labels = {
    'energy_min_steep.xvg': ('Steepest Descent Minimization', 'timesteps', 'energy (kJ / mol)'),
    'energy_min_cg.xvg': ('Conjugate Gradient Minimization', 'timesteps', 'energy (kJ / mol)'),
    'temperature.xvg': ('Temperature', 'time (ps)', 'temperature (K)'),
    'density.xvg': ('Density', 'time (ps)', 'density (g / cm^3)'),
    'rmsd_xray.xvg': ('Root Mean Square w.r.t. Xray', 'time (ps)', 'rmsd (nm)'),
    'rmsd_equil.xvg': ('Root Mean Square w.r.t. Equil', 'time (ps)', 'rmsd (nm)'),
    'radius_gyration.xvg': ('Radius of Gyration', 'time (ps)', 'radius (nm)'),
}


def update_plots(fig: matplotlib.pyplot.Figure, axes2d: List[List[matplotlib.pyplot.Axes]]) -> None:
    """Update the previously initialized plots (if there is any new data).

    Args:
        fig (matplotlib.pyplot.Figure): The main figure
        axes2d (List[List[matplotlib.pyplot.Axes]]): The subplots axes
    """
    nrows = len(axes2d)
    ncols = len(axes2d[0])
    # Declare global variables locally
    global data_glob
    global data_glob_changed
    if data_glob_changed:
        #print('updating plots')
        idx_ax = 0
        for path, data in data_glob:
            #print(num y cols, len(data[0]))
            # Plot multiple y values individually
            for idx_y_col in range(1, len(data[0])):
                if len(data.shape) == 2:
                    xs = data[:,0]
                    ys = data[:,idx_y_col] 
                elif len(data.shape) == 1:
                    #print('path', path)
                    #print('data', data)
                    ys = data[:]
                    xs = np.array([x for x in range(len(ys))])
                else:
                    print('Error! data.shape is not of length 1 or 2:', data.shape)
                    print('path', path)
                    continue
                #print(path)
                #print(xs)
                #print(ys)

                plot_histograms_bool = 'rmsd' in Path(path).stem
                plot_histograms_list = [False, True] if plot_histograms_bool else [False]

                for plot_histogram in plot_histograms_list:
                    # Figure out the axes
                    idx_ax_col = int(idx_ax / ncols)
                    idx_ax_row = idx_ax - (ncols * idx_ax_col)
                    ax = axes2d[idx_ax_col][idx_ax_row]
                    ax.ticklabel_format(style='sci', scilimits=(-2,3), axis='both') # type: ignore
                    if plot_histogram:
                        ax.hist(ys, density=True, histtype='step', bins='sqrt', color='blue') # type: ignore
                        ax.set_xlim(xmin=0)

                        filename = Path(path).stem + '.xvg'
                        if filename in labels:
                            (title, xlabel, ylabel) = labels[filename]
                            ax.set_title(title)
                            ax.set_xlabel(ylabel)
                            ax.set_ylabel('frequency')
                        else:
                            ax.set_title(filename)
                            ax.set_xlabel('y-axis')
                            ax.set_ylabel('frequency')
                    else:
                        ymin = min(ys)
                        ymax = max(ys)
                        ax.set_xlim(min(xs), max(xs))
                        ax.set_ylim(ymin, ymax)

                        # See https://stackoverflow.com/questions/39753282/scatter-plot-with-single-pixel-marker-in-matplotlib
                        ax.scatter(xs, ys, marker='o', s=(72./fig.dpi)**2, c='blue')  # type: ignore

                        # Smooth the timeseries using spline interpolaton
                        spline = UnivariateSpline(xs, ys)
                        spline.set_smoothing_factor(0.15)
                        ys_spline = spline(xs)
                        ax.plot(xs, ys_spline, color='g')

                        # Use Change Point Detection to partition the timeseries
                        # into piecewise 'constant' segments.
                        # Pure python implementation is very slow!
                        # algo = rpt.Pelt(model='rbf', min_size=10).fit(ys)
                        # C implementation is much faster
                        algo = rpt.KernelCPD(kernel='rbf', min_size=10).fit(ys)
                        indices = algo.predict(pen=100) # Large penalty = less segments
                        changepoints = [xs[i-1] for i in indices]
                        for cp_xval in changepoints:
                            ax.vlines(cp_xval, ymin, ymax, color='gray') # type: ignore

                        filename = Path(path).stem + '.xvg'
                        if filename in labels:
                            (title, xlabel, ylabel) = labels[filename]
                            ax.set_title(title)
                            ax.set_xlabel(xlabel)
                            ax.set_ylabel(ylabel)
                        else:
                            ax.set_title(filename)
                            ax.set_xlabel('x-axis')
                            ax.set_ylabel('y-axis')
                    idx_ax += 1
        data_glob_changed = False


def pause_no_show(interval: float) -> None:
    """This is verbatim copy and pasted from matplotlib.pyplot.pause(), with show() commented out.
    show() rudely brings the window to the foreground, which interrupts the user.

    Args:
        interval (float): 
    """
    manager = _pylab_helpers.Gcf.get_active()
    if manager is not None:
        canvas = manager.canvas
        if canvas.figure.stale:
            canvas.draw_idle()
        #show(block=False)
        canvas.start_event_loop(interval)
    else:
        time.sleep(interval)


# NOTE: You should be very careful when using file watchers! Most libraries
# (watchdog, watchfiles, etc) will use operating system / platform-specific
# APIs to check for changes (for performance reasons). However, this can cause
# problems in some cases, specifically for network drives.
# See https://stackoverflow.com/questions/45441623/using-watchdog-of-python-to-monitoring-afp-shared-folder-from-linux
# I'm 99% sure the same problem happens with Docker containers. Either way, the
# solution is to use polling. However, for unknown reasons, simply replacing
# Observer with PollingObserver doesn't seem to work! So we are forced to write
# our own basic file watcher using glob.


"""class TabularDataHandler(PatternMatchingEventHandler):

    def on_any_event(self, event: FileSystemEvent) -> None:
        if event.event_type == 'modified':
            #print(event)
            path = event.src_path
            store_tabular_data(path)"""


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

figure_closed = False

def on_close(event) -> None: # type: ignore
    global figure_closed
    figure_closed = True


def main() -> None:
    global figure_closed
    path = sys.argv[1] if len(sys.argv) > 1 else '.'

    #event_handler = TabularDataHandler(patterns=['*.xvg'])
    #observer = PollingObserver() # This does not work!
    #observer.schedule(event_handler, path, recursive=True)
    #observer.start()

    nrows: int = 3
    ncols: int = 4
    (fig, axes2d) = initialize_plots(nrows, ncols)
    fig.canvas.mpl_connect('close_event', on_close) # type: ignore
    prev_files: Dict[str, float] = {}
    try:
        while not figure_closed:
            # Use our own polling file watcher, see above.
            changed_files = file_watcher_glob(Path(path), '*.xvg', prev_files)
            changed_files_list = list(changed_files.items())
            # Sort by modified time so that if we start RealtimePlots in the
            # middle of the calculation, the plots will still be drawn in the
            # same order.
            changed_files_list.sort(key=lambda x: x[1])
            for file, time_ in changed_files_list:
                print(file)
                store_tabular_data(Path(file))
            prev_files = {**prev_files, **changed_files}

            # NOTE: Do NOT use time.sleep(1.0) here!
            # It does NOT restart the GUI event loop!
            pause_no_show(1.0)  # Wait at least 1 second so we don't just spin.
            update_plots(fig, axes2d)
    except KeyboardInterrupt:
        pass
    #observer.stop()
    #observer.join()

if __name__ == "__main__":
    main()