#!/usr/bin/env python3
import sys
if len(sys.argv) < 3:
    sys.exit("python3 {} <scaf.fa> <N copy [int]>".format(sys.argv[0]))


def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, "".join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, "".join(seq))



def wrapseq(seq):
    l = len(seq)
    for i in range(0,l,80):
        if i+80 <= l:
            print(seq[i:i+80])
        else:
            print(seq[i:])

def main():
    scafs = {}
    need_curated_scafs = []
    with open(sys.argv[1], 'r') as fa:
        for name, seq in read_fasta(fa):
            name = name.replace(">", "").split()[0]
            new_seq = seq * int(sys.argv[2])
        print(">" + name)
        wrapseq(new_seq)


if __name__ == '__main__':
    main()
