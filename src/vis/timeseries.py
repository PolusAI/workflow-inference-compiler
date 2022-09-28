# pylint: disable=C0103, global-statement
import copy
import math
from pathlib import Path
import sys
import time
from typing import List, Dict, Tuple

import matplotlib
from matplotlib import _pylab_helpers # type: ignore
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import UnivariateSpline

import ruptures as rpt

from . import filewatcher

"""from watchdog.observers import Observer
from watchdog.observers.polling import PollingObserver
from watchdog.events import FileSystemEvent, PatternMatchingEventHandler"""


def read_tabular_data(filename: Path) -> List[List[float]]:
    """Reads a tabular data file, skipping comment lines.

    Args:
        filename (Path): The full path of the file to be read

    Returns:
        List[List[float]]: The file contents, with comment lines removed.
    """
    with open(filename, mode='r', encoding='utf-8') as f:
        lines = []
        for line in f.readlines():
            if line.startswith('#'):  # Skip comment lines
                continue
            if line.startswith('@'):  # Skip xmgrace directives
                continue
            lines.append([float(x) for x in line.strip().split()])
    return lines


# Declare global variables to pass data between the main thread and the observer thread. (Not ideal)
data_glob: List[Tuple[Path, np.ndarray]] = []
data_glob_changed: bool = False


def store_tabular_data(filepath: Path, use_stem: bool = True) -> None:
    """Reads the tabular data from filepath and stores it in-memory to be plotted asychronously.

    Args:
        filepath (Path): The tabular data file to be read and stored.
        use_stem (bool, optional): Only store the filename (without extension). Defaults to True.
    """
    # Declare global variables locally
    global data_glob_changed
    floats = read_tabular_data(filepath)

    if len(floats) == 0:
        print('Skipping empty file', filepath)
        return None

    # Check that the array is not ragged; each line must be the same length!
    # I'm not exactly sure why this happens, but it seems like maybe the file
    # contents are not being flushed to disk before getting read back in again.
    # When I manually check the files afterwards, the data is all present.
    lengths = [len(x) for x in floats]
    if not all([length == lengths[0] for length in lengths]):
        print('Warning! Skipping ragged data in', filepath)
        return None

    data = np.array(floats)
    if use_stem:
        filepath = Path(filepath.stem)
    for i , dataglob_i in enumerate(data_glob):
        (pathold, data_old_) = dataglob_i
        if filepath == pathold:
            data_glob[i] = (filepath, data)
            data_glob_changed = True
            return None
    data_glob.append((filepath, data))
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
    # The purpose of calling suptitle is to take up some blank space. suptitle
    # displays a title above all of the subplots, but NOT for the window itself.
    fig.suptitle('')
    plt.get_current_fig_manager().set_window_title('Real-time Analysis Plots') # type: ignore
    fig.tight_layout()
    plt.subplots_adjust(left=0.05, wspace=0.24, bottom=0.05, hspace=0.24)
    plt.show(block=False) # type: ignore
    return (fig, axes2d)


# For now, hardcode labels for the gromacs example.
labels = {
    'energy_min_steep.xvg': ('Steepest Descent Minimization', 'timesteps', 'energy (kJ / mol)'),
    'energy_min_cg.xvg': ('Conjugate Gradient Minimization', 'timesteps', 'energy (kJ / mol)'),
    'temperature.xvg': ('Temperature', 'time (ps)', 'temperature (K)'),
    'density.xvg': ('Density', 'time (ps)', 'density (kg / m^3)'),
    'density.dat': ('Density', 'time (ps)', 'density (g / cm^3)'),
    'energy_total.xvg': ('Total Energy', 'time (ps)', 'energy (kJ / mol)'),
    'rmsd_xray_mainchain.xvg': ('Mainchain RMSD w.r.t. Xray', 'time (ps)', 'rmsd (nm)'),
    'rmsd_equil_mainchain.xvg': ('Mainchain RMSD w.r.t. Equil', 'time (ps)', 'rmsd (nm)'),
    'rmsd_equil_sidechain.xvg': ('Sidechain RMSD w.r.t. Equil', 'time (ps)', 'rmsd (nm)'),
    'radius_gyration.xvg': ('Radius of Gyration', 'time (ps)', 'radius (nm)'),
    'rmsd_equil_ligand_fit.xvg': ('Ligand RMSD (fit) w.r.t. Equil', 'time (ps)', 'rmsd (nm)'),
    'rmsd_equil_ligand_no_fit.xvg': ('Ligand RMSD (no fit) w.r.t. Equil', 'time (ps)', 'rmsd (nm)'),
}


def zscore(mean1: float, std1: float, mean2: float, std2: float) -> float:
    """See https://en.wikipedia.org/wiki/Standard_score"""
    return abs(mean1 - mean2) / math.sqrt(std1**2 + std2**2)


def cluster_intervals_zscore(intervals_0: List[Tuple[int, int]], ys: np.ndarray,
                             zscore_cutoff: float = 1) -> List[List[Tuple[int, int]]]:
    """Performs additional zscore-distance-based clustering of timeseries data
    which has already been initially segmented using Change Point Detection
    (via the `ruptures` library). The purpose of this extra step is to reduce
    the dependence of the clusters on the segmentation penalty threshold parameter.

    Args:
        intervals_ (List[Tuple[int, int]]): The indices which define the intervals
        in which ys has already been segmented.
        ys (np.ndarray): The original timeseries data.
        zscore_cutoff (float, optional): The test statistic threshold used to
        determine additional clustering. Defaults to 1.

    Returns:
        List[List[Tuple[int, int]]]: An updated list of intervals whose
        corresponding y-data are clustered w.r.t. the zscore distance.
    """
    intervals: List[Tuple[int, int]] = copy.deepcopy(intervals_0)
    intervals.sort(key=lambda x: x[1] - x[0], reverse=True) # length of interval

    ys_segmented: Dict[Tuple[int, int], np.ndarray] = {}
    for (i1, i2) in intervals:
        ys_segmented[(i1, i2)] = ys[i1:i2]

    # Greedily cluster intervals based on Z-score
    # It hopefully shouldn't matter, but start with the biggest intervals first
    # (since more I.I.D. data should have better-estimated standard deviations).
    clusters_indices: List[List[Tuple[int, int]]] = [[intervals[0]]]
    init_cluster = ys_segmented[intervals[0]]
    clusters_stats: List[Tuple[float, float]] = [(np.mean(init_cluster), np.std(init_cluster))]
    for i1, i2 in intervals[1:]:
        ys_seg = ys_segmented[(i1, i2)]

        z_scores = []
        for i, (mean_c, std_c) in enumerate(clusters_stats):
            z_score = zscore(mean_c, std_c, np.mean(ys_seg), np.std(ys_seg))
            z_scores.append((i, z_score))
        z_scores.sort(key=lambda x: x[1])
        zscore_seg = z_scores[0][1]
        if zscore_seg < zscore_cutoff:
            # Statistically equivalent, so add segment to existing cluster.
            cluster_index = z_scores[0][0]
            clusters_indices[cluster_index].append((i1, i2))
        else:
            # Create a new cluster
            clusters_indices.append([(i1, i2)])

        # Update stats after each iteration
        clusters_stats = []
        for cluster_indices in clusters_indices:
            ys_cluster = [ys_segmented[i1_i2] for i1_i2 in cluster_indices]
            ys_cluster_flat = np.concatenate(ys_cluster)
            cluster_stats = (np.mean(ys_cluster_flat), np.std(ys_cluster_flat))
            clusters_stats.append(cluster_stats)

    return clusters_indices


def update_plots(fig: matplotlib.pyplot.Figure, axes2d: List[List[matplotlib.pyplot.Axes]]) -> None:
    """Update the previously initialized plots (if there is any new data).

    Args:
        fig (matplotlib.pyplot.Figure): The main figure
        axes2d (List[List[matplotlib.pyplot.Axes]]): The subplots axes
    """
    nrows = len(axes2d)
    ncols = len(axes2d[0])
    # Declare global variables locally
    global data_glob_changed
    if data_glob_changed:
        #print('updating plots')
        idx_ax = 0
        for path, data in data_glob:
            if idx_ax >= nrows*ncols:
                # TODO: resize plot grid
                # initialize_plots(...)
                break
            #print(num y cols, len(data[0]))
            # Plot multiple y values individually
            col_indices = list(range(1, len(data[0])))
            if 'radius_gyration' == Path(path).stem:
                col_indices = [1]
            for idx_y_col in col_indices:
                if len(data.shape) == 2:
                    xs = data[:,0]
                    ys = data[:,idx_y_col]
                elif len(data.shape) == 1:
                    #print('path', path)
                    #print('data', data)
                    ys = data[:]
                    xs = np.array(list(range(len(ys))))
                else:
                    print('Error! data.shape is not of length 1 or 2:', data.shape)
                    print('path', path)
                    continue
                #print(path)
                #print(xs)
                #print(ys)

                interval_indices_clustered = [[(0, len(xs))]]
                changepoints = []
                min_size=10
                if len(xs) > min_size:
                    # Use Change Point Detection to partition the timeseries
                    # into piecewise 'constant' segments.
                    # Pure python implementation is very slow!
                    # algo = rpt.Pelt(model='rbf', min_size=10).fit(ys)
                    # C implementation is much faster
                    try:
                        algo = rpt.KernelCPD(kernel='rbf', min_size=min_size).fit(ys)
                        indices = algo.predict(pen=100) # Large penalty = less segments
                        indices_zero = [0] + indices
                        interval_indices = list(zip(indices_zero, indices_zero[1:]))
                        interval_indices_clustered = cluster_intervals_zscore(interval_indices, ys)
                        changepoints = [xs[i-1] for i in indices]
                    except rpt.exceptions.BadSegmentationParameters as bsp:
                        pass

                xs_segmented = []
                ys_segmented = []
                for inter in interval_indices_clustered:
                    xs_seg_nest = [xs[i1:i2] for i1, i2 in inter]
                    xs_seg = [x for y in xs_seg_nest for x in y]
                    xs_segmented.append(xs_seg)
                    ys_seg_nest = [ys[i1:i2] for i1, i2 in inter]
                    ys_seg = [x for y in ys_seg_nest for x in y]
                    ys_segmented.append(ys_seg)

                plot_histograms_bool = 'rmsd' in Path(path).stem
                plot_histograms_list = [False, True] if plot_histograms_bool else [False]

                for plot_histogram in plot_histograms_list:
                    # Figure out the axes
                    idx_ax_col = int(idx_ax / ncols)
                    idx_ax_row = idx_ax - (ncols * idx_ax_col)
                    ax = axes2d[idx_ax_col][idx_ax_row]
                    ax.clear() # type: ignore
                    ax.ticklabel_format(style='sci', scilimits=(-2,3), axis='both') # type: ignore
                    cmap = plt.get_cmap("tab10") # type: ignore
                    colors = [cmap(i) for i in range(len(xs_segmented))]

                    # For plotting purposes only, ignore first 1% of data
                    # (Initial data may have transients that distort plots)
                    ymin = min(ys[int(len(ys)/100):])
                    ymax = max(ys[int(len(ys)/100):])
                    if plot_histogram:
                        if not (np.isnan(ymin) or np.isnan(ymax)):
                            ax.hist(ys_segmented, range=(ymin, ymax), stacked=True, density=True, # type: ignore
                                    histtype='barstacked', bins='sqrt', color=colors) # type: ignore
                        else:
                            ax.hist(ys_segmented, stacked=True, density=True, # type: ignore
                                    histtype='barstacked', bins='sqrt', color=colors) # type: ignore
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
                        ax.set_xlim(min(xs), max(xs))
                        if not (np.isnan(ymin) or np.isnan(ymax)):
                            ax.set_ylim(ymin, ymax)

                        # See https://stackoverflow.com/questions/39753282/scatter-plot-with-single-pixel-marker-in-matplotlib
                        for i, (xs_seg, ys_seg) in enumerate(zip(xs_segmented, ys_segmented)):
                            ax.scatter(xs_seg, ys_seg, marker='o', s=(72./fig.dpi)**2, color=colors[i])  # type: ignore

                        # Smooth the timeseries using spline interpolaton of degree k
                        k=3
                        if len(ys) > k:
                            spline = UnivariateSpline(xs, ys, k=k)
                            # "Number of knots will be increased until the smoothing condition is satisfied:
                            # sum((w[i] * (y[i]-spl(x[i])))**2, axis=0) <= s"
                            # The following is somewhat ad-hoc (based on the above condition)
                            # but currently produces visually appealing results.
                            std_ys = np.std(ys)
                            sfactor = abs(np.mean(ys))
                            if std_ys < 1.0: # if 'rmsd' in Path(path).stem or 'gyration' in Path(path).stem:
                                sfactor = math.sqrt(sfactor)
                            if 'total' in Path(path).stem:
                                # The total energy is controlled by the thermostat,
                                # so its std will be relatively smaller w.r.t.
                                # the other plots, so square it to make it bigger.
                                sfactor = (std_ys ** 2) * sfactor
                            else:
                                sfactor = math.sqrt(std_ys) * sfactor
                            if std_ys < 1.0: # if 'rmsd' in Path(path).stem or 'gyration' in Path(path).stem:
                                sfactor = math.sqrt(sfactor)
                            #print('path, sfactor', path, sfactor)
                            spline.set_smoothing_factor(sfactor)
                            ys_spline = spline(xs)
                            ax.plot(xs, ys_spline, color='white')

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
        interval (float): Delay execution for a given number of seconds.
        The argument may be a floating point number for subsecond precision.
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


"""class TabularDataHandler(PatternMatchingEventHandler):

    def on_any_event(self, event: FileSystemEvent) -> None:
        if event.event_type == 'modified':
            #print(event)
            path = event.src_path
            store_tabular_data(path)"""


figure_closed = False

def on_close(event) -> None: # type: ignore
    """This event handler is used to terminate the polling loop when the matplotlib figure is closed."""
    global figure_closed
    figure_closed = True


def main() -> None:
    """See udocs/userguide.md#real-time-plots"""
    cachedir_path = sys.argv[1] if len(sys.argv) > 1 else '.'
    file_patterns = sys.argv[2:] if len(sys.argv) > 2 else ['*.xvg', '*.dat']

    #event_handler = TabularDataHandler(patterns=[file_pattern])
    #observer = PollingObserver() # This does not work!
    #observer.schedule(event_handler, cachedir_path, recursive=True)
    #observer.start()

    nrows: int = 3
    ncols: int = 4
    (fig, axes2d) = initialize_plots(nrows, ncols)
    fig.canvas.mpl_connect('close_event', on_close) # type: ignore
    prev_files: Dict[str, float] = {}
    try:
        while not figure_closed:
            # Use our own polling file watcher, see above.
            changed_files = filewatcher.file_watcher_glob(Path(cachedir_path), file_patterns, prev_files)
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
