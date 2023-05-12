#!/usr/bin/env python3
import sys
import os

if len(sys.argv) < 3:
    print("Usage: python3 {} <asm2ref.paf> <paf|fmash>\n".format(sys.argv[0]))
    exit()
#elif sys.argv[2] != "paf" and sys.argv[2] != "fmash":
#    print("only recognize paf or fmash (mashmap like) file format")

records = {}
suffix = sys.argv[2]

with open(sys.argv[1], 'r') as fh:
    for i in fh:
        tmp = i.strip().split()
        ref = tmp[5]
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
    with open(outdir + "/" + k + "." + suffix, 'w') as out:
        for i in records[k]:
            out.write(" ".join(i))
            out.write("\n")
