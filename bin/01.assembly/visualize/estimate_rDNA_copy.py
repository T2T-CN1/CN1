# !/usr/bin/env python3
# -*- coding: utf-8 -*-

import seaborn as sns
from icecream import ic
import numpy as np
import pandas as pd
import sys

from matplotlib import (pyplot as plt, lines)

def parse_depth(depth_input, coverage):
    kmer_copies = []
    total_count = []
    with open(depth_input) as fh:
        for row in fh:
            pos, kmerseq, count = row.strip().split()
            kmer_copies.append([int(pos), round(int(count) / coverage)],)
            total_count.append(round(int(count) / coverage))
    average = sum(total_count) / len(kmer_copies)
    df = pd.DataFrame(data=kmer_copies,
                 columns = ['pos','copies'],
                 index = range(0, len(kmer_copies)))

    return df, average


def plot_depth(df, average, output, max_length=None, max_copy=None):
    plot_title = "estimated rDNA copies"
    plt.title(plot_title)
    fig, ax = plt.subplots(figsize = (15, 4))
    g = sns.lineplot(data=df, x="pos", y="copies")
    g.set(xlabel='Genome Position (bp)', ylabel="Copies", ylim=(0, max_copy), xlim=(0, max_length))

    # add average of depth line
    # gene features
    features = {(3658,5526): '18S',
                (6597,6753): '5.8S',
                (7921,12971): '28S'}
    for region in features:
        start, end = region
        gene_name = features[region]
        ax.add_line(lines.Line2D([start], [df['copies'][start-1], 0], color='purple', ls='-', lw=0.2))
        ax.add_line(lines.Line2D([end], [df['copies'][end-1], 0], color='purple', ls='-', lw=0.2))
        x1, y1 =  (start + end) /2 , average * 1.2
        ax.text(x1, y1, gene_name, ha="center")


    #plt.figure(figsize=(8, 4))
    # bbox_inches='tight'
    plt.savefig(output, format='pdf')
    plt.close()


if __name__ == '__main__':
    if len(sys.argv) < 4:
        sys.exit(f"python3 {sys.argv[0]} rdnaRef.31mer.freq_in_rawreads.txt output coverage")
    else:
        coverage = int(sys.argv[3])
        output = sys.argv[2]
        kmer_copies, average = parse_depth(sys.argv[1], coverage)
        plot_depth(kmer_copies, average, output, max_length=14000, max_copy=300)
