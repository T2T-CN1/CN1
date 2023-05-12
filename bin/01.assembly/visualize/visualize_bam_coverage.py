# !/usr/bin/env python3
# -*- coding: utf-8 -*-

import seaborn as sns
import os
import numpy as np
from icecream import ic
import argparse

from matplotlib import (pyplot as plt, lines)

def parse_depth(depth_input):
    """Parse depth file.

    Args:
        depth_input (str): Path to depth file.
        genome_size (int): Genome size.

    Returns:
        list: List with depth.

    """
    #length = os.system('wc -l ' + depth_input + " |awk '{print $1}' ")
    #ic(length)
    #print(type(length))
    # no more than 250 Mb
    length = 250000000
    depth = [0] * int(length)
    values = []
    pos = []
    references = set()
    with open(depth_input) as depth_object:
        for row in depth_object:
            genome_id, position, depth_val = row.split()

            references.add(genome_id)

            if len(references) > 1:
                raise Exception(' This script only handles one genome - contig.')

            if float(depth_val) > args.maxdepth:
                depth_val = args.maxdepth

            depth[int(position)] = float(depth_val)
            values.append(float(depth_val))
            pos.append(int(position))
    return depth, pos, values


def plot_depth(depth_report, annotation, output_name, plot_title, normalize=False, depth_cut_off=5):

    data, pos, values = parse_depth(depth_report)
    min_position = min(pos)
    max_position = max(pos)

    if normalize:
        data = [xx / max(data) for xx in data]
        y_label = "Normalized Depth"
    else:
        y_label = "Depth"
    sns.set(color_codes=True)
    plt.title(plot_title)
    #fig, ax = plt.subplot(111)
    fig, ax = plt.subplots(figsize = (8, 4))

    sns_plot = sns.lineplot(x=pos, y=values, lw=0.3, color='grey')
    sns_plot.set(xlabel='Genome Position (bp)', ylabel=y_label)

    if not normalize:
        ax.add_line(lines.Line2D([min_position, max_position + 1], [depth_cut_off], color="r", ls='--', dash_capstyle='butt', lw=0.3))

    # Annotation regions
    colors = ['blue', 'green', 'red', 'orange', 'yellow']
    if annotation:
        index = 0
        with open(annotation, 'r') as fh:
            for bed in fh:
                depth_this_region = []
                scaf, start, end, annot = bed.strip().split("\t")
                start = int(start)
                end = int(end)
                middle_pos = round((start + end) / 2)
                for p in range(start, end):
                    ic(p, data[p])
                    depth_this_region.append(data[p])
                average_depth = sum(depth_this_region) / len(depth_this_region)
                ic(average_depth)
                # add average of depth line
                ax.add_line(lines.Line2D([start, end], [average_depth], color=colors[index], ls='-', lw=0.8))
                x1, y1 =  middle_pos, average_depth * 1.5
                text = annot
                ax.text(x1, y1, text, ha="center")
                index += 1


    #plt.figure(figsize=(8, 4))
    # bbox_inches='tight'
    plt.savefig(output_name, format='pdf')
    plt.close()

    print("Done :)")


if __name__ == '__main__':
    usage = """Plot genome Depth across genome.

    Args:
        depth_report (str): Path to samtool's depth file.
        output_name (str): Path to output PNG image.
        plot_title (str): Plot title.
        genome_size (int): Genome size.
        normalize (bool): If `True`, normalizes the depth by the largest depth (default = `False`).
        depth_cut_off (int): Plot a line to represent a targeted depth (default = 20).

    """
    parser = argparse.ArgumentParser(
    description=usage,
    formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('depth', metavar='<File>', type=str,
                        help="input bam depth")

    parser.add_argument('-b', metavar='<BED>', type=str, required=False,
                        dest='annotation', help="different regions")

    parser.add_argument('-o', metavar='<STR>', type=str, required=True,
                        dest='output', help="output name")

    parser.add_argument('-t', metavar='<STR>', type=str, required=False,
                        dest='title', help="title for plotting")

    parser.add_argument('-n', action="store_true",
                        dest='normalize', help="whether to normalize depth")

    parser.add_argument('-d', metavar='<INT>', type=int,
                        dest='depth_cutoff', help="depth cutoff for plotting default=5",
                        default=5)

    parser.add_argument('-mp', metavar='<INT>', type=int,
                        dest='maxdepth', help="maximum depth for plotting default=300",
                        default=300)

    args = parser.parse_args()
    plot_depth(args.depth, args.annotation, args.output, args.title, normalize=args.normalize, depth_cut_off=args.depth_cutoff)
