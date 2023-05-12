#!/usr/bin/env python3

import sys
from icecream import ic
import os

if len(sys.argv) < 2:
    sys.exit(f"python3 {sys.argv[0]} population matrix")

def main():
    sample2popu = {}
    sample_order = {}
    with open(sys.argv[1], 'r') as fh:
        for i in fh:
            tmp = i.strip().split()
            sample, subpopu, popu = tmp
            sample2popu[sample] = popu


    with open(sys.argv[2], 'r') as fh:
        for i in fh:
            tmp = i.strip().split()
            if tmp[0].startswith("#"):
                for index, sample in enumerate(tmp[5:]):
                    feild = index + 5
                    sample = sample.split(".")[0]
                    sample_order[feild] = sample
            else:
                popu_score = {}
                pid = "\t".join(tmp[0:5])
                for index, info in enumerate(tmp[5:]):
                    feild = index + 5
                    sample = sample_order[feild]
                    popu = sample2popu[sample]
                    if info == ".":
                        if popu not in popu_score:
                            popu_score[popu] = [0,]
                        else:
                            popu_score[popu].append(0)
                    else:
                        if popu not in popu_score:
                            popu_score[popu] = [1,]
                        else:
                            popu_score[popu].append(1)

                output = [pid,]
                for p in sorted(popu_score.keys()):
                    p_len = len(popu_score[p])
                    exists= popu_score[p].count(1)
                    not_exists = popu_score[p].count(0)
                    output.append(f"{p}:{p_len}:{exists}:{not_exists}")
                print("\t".join(output))


if __name__ == '__main__':
    main()
