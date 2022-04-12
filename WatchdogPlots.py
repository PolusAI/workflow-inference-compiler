from pathlib import Path
import sys
import time
from typing import Any

import numpy as np
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


# Create blank placeholder subplots
nrows: int = 2
ncols: int = 3
plt.style.use('dark_background')
fig, axes2d = plt.subplots(nrows=nrows, ncols=ncols, figsize=(12, 9))
fig.suptitle('')
fig.tight_layout()
plt.subplots_adjust(left=0.08, wspace=0.24, bottom=0.08, hspace=0.24)
# NOTE: "The pause is needed because the GUI events happen while the main code is sleeping, including drawing."
# See https://stackoverflow.com/questions/28269157/plotting-in-a-non-blocking-way-with-matplotlib
plt.pause(0.1)
plt.show(block=False) # type: ignore


# Declare global variables to pass data between the main thread and the observer thread. (Not ideal)
idx_glob = -1
data_glob = None
path_glob = None


def plot_new_data():
    # Declare global variables locally
    global idx_glob
    global data_glob
    global path_glob
    if isinstance(data_glob, np.ndarray):
        data = data_glob
        data_glob = None
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
                    ax.ticklabel_format(style='sci', scilimits=(-2,3), axis='both')
                    # See https://stackoverflow.com/questions/39753282/scatter-plot-with-single-pixel-marker-in-matplotlib
                    ax.scatter(xs, ys, marker='o', s=(72./fig.dpi)**2)  # type: ignore
                    ax.set_title(Path(path_glob).stem + '.xvg')
                    ax.set_xlabel('x-axis')
                    ax.set_ylabel('y-axis')
                    plt.draw()
                    plt.pause(0.1)
                idx += 1


class XvgHandler(PatternMatchingEventHandler):

    def __init__(self):
        patterns = ['*.xvg']
        super().__init__(patterns)

    def on_any_event(self, event: FileSystemEvent) -> Any:
        # Declare global variables locally
        global idx_glob
        global data_glob
        global path_glob
        if event.event_type == 'modified':
            #print(event)
            idx_glob += 1
            data_glob = read_tabular_data(event.src_path)
            path_glob = event.src_path


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else '.'
    event_handler = XvgHandler()
    observer = Observer()
    observer.schedule(event_handler, path, recursive=True)
    observer.start()
    try:
        while True:
            time.sleep(1)
            plot_new_data()
    except KeyboardInterrupt:
        observer.stop()
    observer.join()


if __name__ == "__main__":
    main()