#!/usr/bin/env python3
import sys
if len(sys.argv) < 3:
    sys.exit(f"python3 {sys.argv[0]} allscaf.lens linkview.order.txt")


lens = {}
with open(sys.argv[1], 'r') as fh:
    for i in fh:
        tmp = i.strip().split()
        lens[tmp[0]] = tmp[1]

with open(sys.argv[2], 'r') as fh:
    for i in fh:
        tmp = i.strip().split()
        for index,scaf in enumerate(tmp):
            if scaf in lens:
                tmp[index] = f"{scaf}:0:{lens[scaf]}"
        print(" ".join(tmp))


