#!/usr/bin/env python3
import sys
import os


records = {}

with open(sys.argv[1], 'r') as fh:
    for i in fh:
        tmp = i.strip().split()[:3]
        ref = tmp[0]
        if ref not in records.keys():
            records[ref] = [tmp,]
        else:
            records[ref].append(tmp)
if not os.path.exists("splitbyChrs"):
    os.mkdir("splitbyChrs")


for k in sorted(records.keys()):
    outdir = "splitbyChrs/" + k
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    with open(outdir + "/rdna.bed.txt", 'w') as out:
        for i in records[k]:
            out.write(" ".join(i) + " red")
            out.write("\n")
