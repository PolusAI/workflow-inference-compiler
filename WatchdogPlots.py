from pathlib import Path
import sys
import time
from typing import List, Tuple
import matplotlib

import numpy as np
from matplotlib import _pylab_helpers # type: ignore
import matplotlib.pyplot as plt

from watchdog.observers import Observer
from watchdog.events import FileSystemEvent, PatternMatchingEventHandler


def read_tabular_data(filename: str) -> np.ndarray:
    lines = []
    for line in open(filename, 'r').readlines():
        if line.startswith('#'):  # Skip comment lines
            continue
        lines.append([float(x) for x in line.strip().split()])
    arr = np.array(lines)
    return arr


# Declare global variables to pass data between the main thread and the observer thread. (Not ideal)
idx_glob = -1
data_glob: np.ndarray = np.array([])
path_glob: Path


def initialize_plots(nrows: int, ncols: int) -> Tuple[matplotlib.pyplot.Figure, List[List[matplotlib.pyplot.Axes]]]:
    """Initialize a grid of blank subplots.

    Args:
        nrows (int): the number of rows
        ncols (int): the number of columns

    Returns:
        Tuple[matplotlib.pyplot.Figure, List[List[matplotlib.pyplot.Axes]]]: The main figure and subplots axes
    """
    plt.style.use('dark_background') # type: ignore
    fig, axes2d = plt.subplots(nrows=nrows, ncols=ncols, figsize=(12, 9))
    fig.suptitle('')
    fig.tight_layout()
    plt.subplots_adjust(left=0.08, wspace=0.24, bottom=0.08, hspace=0.24)
    plt.show(block=False) # type: ignore
    return (fig, axes2d)


def update_plots(fig: matplotlib.pyplot.Figure, axes2d: List[List[matplotlib.pyplot.Axes]]) -> None:
    """Update the previously initialized plots (if there is any new data).

    Args:
        fig (matplotlib.pyplot.Figure): The main figure
        axes2d (List[List[matplotlib.pyplot.Axes]]): The subplots axes
    """
    nrows = len(axes2d)
    ncols = len(axes2d[0])
    # Declare global variables locally
    global idx_glob
    global data_glob
    global path_glob
    if data_glob.size > 0:
        data = data_glob
        data_glob = np.array([])
        idx = 0
        for axes1d in axes2d:
            for ax_ in axes1d:
                if idx == idx_glob:
                    idx_col = int(idx_glob / ncols)
                    idx_row = idx_glob - (ncols * idx_col)
                    ax = axes2d[idx_col][idx_row]
                    xs = data[:,0]
                    ys = data[:,1]
                    #print(path_glob)
                    #print(xs)
                    #print(ys)
                    ax.set_xlim(min(xs), max(xs))
                    ax.set_ylim(min(ys), max(ys))
                    ax.ticklabel_format(style='sci', scilimits=(-2,3), axis='both') # type: ignore
                    # See https://stackoverflow.com/questions/39753282/scatter-plot-with-single-pixel-marker-in-matplotlib
                    ax.scatter(xs, ys, marker='o', s=(72./fig.dpi)**2)  # type: ignore
                    ax.set_title(path_glob.stem + '.xvg')
                    ax.set_xlabel('x-axis')
                    ax.set_ylabel('y-axis')
                idx += 1


class XvgHandler(PatternMatchingEventHandler):

    def __init__(self) -> None:
        patterns = ['*.xvg']
        super().__init__(patterns)

    def on_any_event(self, event: FileSystemEvent) -> None:
        # Declare global variables locally
        global idx_glob
        global data_glob
        global path_glob
        if event.event_type == 'modified':
            #print(event)
            idx_glob += 1
            data_glob = read_tabular_data(event.src_path)
            path_glob = Path(event.src_path)


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


def main() -> None:
    path = sys.argv[1] if len(sys.argv) > 1 else '.'
    event_handler = XvgHandler()
    observer = Observer()
    observer.schedule(event_handler, path, recursive=True)
    observer.start()

    nrows: int = 2
    ncols: int = 3
    (fig, axes2d) = initialize_plots(nrows, ncols)
    try:
        while True:
            # NOTE: Do NOT use time.sleep(0.1) here!
            # It does NOT restart the GUI event loop!
            pause_no_show(0.1)
            update_plots(fig, axes2d)
    except KeyboardInterrupt:
        observer.stop()
    observer.join()


if __name__ == "__main__":
    main()