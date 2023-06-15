#!/usr/bin/env python3

import sys
import os
from icecream import ic

if len(sys.argv) < 4:
    sys.exit(f"python3 {sys.argv[0]} *.bed density_cutoff outprefix")


def rank_heterozygousity(val):
    if val > 0 and val < 2:
        return 1
    elif val >= 2 and val < 5:
        return 2
    elif val >= 5 and val < 10:
        return 3
    elif val >= 10 and val < 15:
        return 4
    elif val >= 15 and val < 20:
        return 5
    elif val >= 20:
        return 6

#def rank_heterozygousity(value, threshold):
    #step = 0.00
    #if threshold >= step:
     #   sys.exit(f"the threshold must less than step (0.002)")
    #if value > threshold and value < 0.002:
    #    rank = 1
    #else:
    #    rank = value // step + 1
    #if rank > 6:
    #    rank = 6
    #rank = str(rank).split(".")[0]
    #return rank



def add_segment(all_seg, blocks, status, threshold):
    heterozygousity = []
    for l in blocks:
        heterozygousity.append(float(l.split()[4]))
    average_heterozygousity = sum(heterozygousity) / len(heterozygousity)
    rank = rank_heterozygousity(average_heterozygousity)
    chrom = blocks[0].split()[0]
    start = blocks[0].split()[1]
    end = blocks[-1].split()[2]
    segment_id = f"{chrom}:{start}:{end}"
    if chrom not in all_seg:
        all_seg[chrom] = [(segment_id, status, rank), ]
    else:
        all_seg[chrom].append((segment_id, status, rank))

def get_status(density, threshold):
    if density >= threshold:
        return "hete"
    else:
        return "homo"

def main():
    threshold = float(sys.argv[2])
    outpre = sys.argv[3]
    block_sequence = {} # key is chrom, val is block sequences
    link = set()

    """
    chr1    0   1000000 10  10  1000000 0.0000100
    chr1    1000000 2000000 11  11  1000000 0.0000110
    """

    chrom = "" # condition 1 of stopping chain
    # hete/homo/ condition 2 of stopping chain
    blocks = [] # tmp collection
    with open(sys.argv[1], 'r') as fh:
        for line in fh:
            # the first line
            tmp = line.strip().split("\t")
            density = int(tmp[4]) # count
            if chrom == "":
                blocks = [line,]
                chrom = tmp[0]
                current_status = get_status(density, threshold)
            # still in same chrom
            elif tmp[0] == chrom:
                if get_status(density, threshold) == current_status:
                    blocks.append(line)
                else: # hete/home type is changed
                    add_segment(block_sequence, blocks, current_status, threshold)
                    blocks = [line,]
                    current_status = get_status(density, threshold)
            else: # next chrom
                chrom = tmp[0]
                add_segment(block_sequence, blocks, current_status, threshold)
                blocks = [line,]
                current_status = get_status(density, threshold)

        add_segment(block_sequence, blocks, current_status, threshold)  # the last one

    out = open(outpre + ".node.color.csv", 'w')
    gfa = open(outpre + ".graph.gfa", 'w')
    for chrom in block_sequence:
        this_chr_block_sequence = block_sequence[chrom]
        links = []
        merged_block_len = len(this_chr_block_sequence)
        tagged_segment = set()
        for i in range(0, merged_block_len-1):
            j = i + 1
            left_sid, left_tag, left_hete_rank = this_chr_block_sequence[i]
            right_sid, right_tag, right_hete_rank = this_chr_block_sequence[j]
            # only two conditions, homo -> hete, hete -> home
            if left_tag == "homo":
                l1 = f"L\tc{i}_{left_sid}\t+\t{right_hete_rank}m{j}_{right_sid}\t+\t100M"
                l2 = f"L\tc{i}_{left_sid}\t+\t{right_hete_rank}p{j}_{right_sid}\t+\t100M"
                links.append(l1)
                links.append(l2)
                tagged_segment.add(f"c{i}_{left_sid}")
                tagged_segment.add(f"{right_hete_rank}m{j}_{right_sid}")
                tagged_segment.add(f"{right_hete_rank}p{j}_{right_sid}")
            else:
                l1 = f"L\t{left_hete_rank}m{i}_{left_sid}\t+\tc{j}_{right_sid}\t+\t100M"
                l2 = f"L\t{left_hete_rank}p{i}_{left_sid}\t+\tc{j}_{right_sid}\t+\t100M"
                links.append(l1)
                links.append(l2)
                tagged_segment.add(f"{left_hete_rank}m{i}_{left_sid}")
                tagged_segment.add(f"{left_hete_rank}p{i}_{left_sid}")
                tagged_segment.add(f"c{j}_{right_sid}")


        # output segment
        # color label
        rank2color = {
            '1p':'#f7fbff',
            '2p':'#deebf7',
            '3p':'#c6dbef',
            '4p':'#9ecae1',
            '5p':'#6baed6',
            '6p':'#4292c6',
            '1m':'#fff5f0',
            '2m':'#fee0d2',
            '3m':'#fcbba1',
            '4m':'#fc9272',
            '5m':'#fb6a4a',
            '6m':'#ef3b2c',

        }
        for s in tagged_segment:
            tmp = s.split(":")
            length = int(tmp[2]) - int(tmp[1])
            node = tmp[0]
            if node.startswith("c"):
                color = "#B8B8B4"
            else:
                color = rank2color[node[:2]]
            print(f"{s},{color}", file=out)
            print(f"S\t{s}\t*\tLN:i:{length}", file=gfa)
        # output links
        for l in links:
            print(l, file=gfa)
    out.close()
    gfa.close()


if __name__ == '__main__':
    main()
