#!/usr/bin/env python3
import sys
if len(sys.argv) < 3:
    sys.exit("python3 {} <scaf.fa> <*.strand>".format(sys.argv[0]))


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


def revcom_transtable(sequence):
    # make a sequence complement #
    transtable = str.maketrans('ATCGatcg', 'TAGCtagc')
    sequence = sequence.translate(transtable)
    return sequence[::-1]


def main():
    with open(sys.argv[2], 'r') as fh:
        ori = {}
        for i in fh:
            [readid, strand] =  i.strip().split()
            ori[readid] = strand
    with open(sys.argv[1], 'r') as fa:
        for name, seq in read_fasta(fa):
            name = name.replace(">", "").split()[0]
            if name in ori and ori[name] == '-':
                print(f">{name}\n{revcom_transtable(seq)}")
            else:
                print(f">{name}\n{seq}")

if __name__ == '__main__':
    main()
