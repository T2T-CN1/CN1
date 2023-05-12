#/usr/bin/env python3
import os
import sys

if len(sys.argv) < 3:
    sys.exit(f"python3 {sys.argv[0]} db.kmer.freq.txt target.kmer.txt")

def main():
    target_pos_kmer = {} # start pos => kmer seq
    target_kmer_freq = {} # kmer seq => count
    # read target list
    with open(sys.argv[2], 'r') as fh:
        for i in fh:
            tmp = i.strip().split()
            pos = tmp[0]
            kmer = tmp[1]
            target_pos_kmer[pos] = kmer
            target_kmer_freq[kmer] = 1
    # scan db
    with open(sys.argv[1], 'r') as fh:
        for i in fh:
            tmp = i.strip().split()
            kmer = tmp[0]
            count = tmp[1]
            if kmer in target_kmer_freq:
                target_kmer_freq[kmer] = count


    for pos in sorted(target_pos_kmer.keys()):
        print(f"{pos}\t{target_pos_kmer[pos]}\t{target_kmer_freq[target_pos_kmer[pos]]}")


if __name__ == '__main__':
    main()
