#!/usr/bin/env python3

import sys
import os
from icecream import ic

if len(sys.argv) < 3:
    sys.exit(f"python3 {sys.argv[0]} *.bed density_cutoff")


def rank_heterozygousity(value, threshold):
    step = 0.002
    if threshold >= step:
        sys.exit(f"the threshold must less than step (0.002)")
    if value > threshold and value < 0.002:
        rank = 1
    else:
        rank = value // step + 1
    if rank > 6:
        rank = 6
    rank = str(rank).split(".")[0]
    return rank



def add_segment(all_seg, blocks, status, threshold):
    heterozygousity = []
    for l in blocks:
        heterozygousity.append(float(l.split()[-1]))
    average_heterozygousity = sum(heterozygousity) / len(heterozygousity)
    rank = rank_heterozygousity(average_heterozygousity, threshold)
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
            density = float(tmp[6])
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

    out = open("rank.txt", 'w')
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
                print(right_hete_rank, file=out)
                l1 = f"L\tc{i}_{left_sid}\t+\t{right_hete_rank}m{j}_{right_sid}\t+\t100M"
                l2 = f"L\tc{i}_{left_sid}\t+\t{right_hete_rank}p{j}_{right_sid}\t+\t100M"
                links.append(l1)
                links.append(l2)
                tagged_segment.add(f"c{i}_{left_sid}")
                tagged_segment.add(f"{right_hete_rank}m{j}_{right_sid}")
                tagged_segment.add(f"{right_hete_rank}p{j}_{right_sid}")
            else:
                print(left_hete_rank, file=out)
                l1 = f"L\t{left_hete_rank}m{i}_{left_sid}\t+\tc{j}_{right_sid}\t+\t100M"
                l2 = f"L\t{left_hete_rank}p{i}_{left_sid}\t+\tc{j}_{right_sid}\t+\t100M"
                links.append(l1)
                links.append(l2)
                tagged_segment.add(f"{left_hete_rank}m{i}_{left_sid}")
                tagged_segment.add(f"{left_hete_rank}p{i}_{left_sid}")
                tagged_segment.add(f"c{j}_{right_sid}")


        # output segment
        for s in tagged_segment:
            tmp = s.split(":")
            length = int(tmp[2]) - int(tmp[1])
            print(f"S\t{s}\t*\tLN:i:{length}")
        # output links
        for l in links:
            print(l)
    out.close()


if __name__ == '__main__':
    main()
