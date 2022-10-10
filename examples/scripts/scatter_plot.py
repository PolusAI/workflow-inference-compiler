import argparse

import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('--output_png_path', default='scatter.png')
parser.add_argument('--xs', nargs='+', default=[], type=float)
parser.add_argument('--ys', nargs='+', default=[], type=float)
parser.add_argument('--ys2', nargs='+', default=[], type=float)
args = parser.parse_args()

plt.scatter(args.xs, args.ys)
if args.ys2 != []:
    plt.scatter(args.xs, args.ys2)

plt.savefig(args.output_png_path)
