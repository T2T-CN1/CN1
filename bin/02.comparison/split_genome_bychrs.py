#!/usr/bin/env python3
import os
import sys
if len(sys.argv) < 3:
    print("Usage: python3 {} {} {}".format(sys.argv[0], "mat.fas", "pat.fas"))
    exit()


def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, '\n'.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, '\n'.join(seq))

for fa in [sys.argv[1], sys.argv[2]]:
    with open(fa, 'r') as fp:
        for name, seq in read_fasta(fp):
            name = name.replace(">", "")
            seq_id = name.split(" ")[0]
            chrid = seq_id.split("_")[1]
            if os.path.exists(chrid) == False:
                os.mkdir(chrid)
            with open(chrid + "/" + seq_id + ".fa", 'w') as fout:
                fout.write(">" + seq_id + "\n" + seq + "\n")
