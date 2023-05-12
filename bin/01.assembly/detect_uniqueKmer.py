#!/usr/bin/env python3 
import os
import sys

if len(sys.argv) < 2:
    sys.exit(f"python3 {sys.argv[0]} a.kmer.txt b.kmer.txt ...")


def main():
    code = {}
    all_kmers = {}
    all_kmers_sum = {}
    for index, kmer_f in enumerate(sys.argv[1:]):
        index = index + 1
        code[index] = kmer_f.replace(".sort.txt", "")
        with open(kmer_f, 'r') as fh:
            for i in fh:
                tmp = i.strip().split()
                kmer_seq = tmp[0]
                if kmer_seq not in all_kmers_sum:
                    all_kmers_sum[kmer_seq] = index
                else:
                    all_kmers_sum[kmer_seq] += index
                all_kmers[kmer_seq] = index

    count = 0
    for k in sorted(all_kmers.keys()):
        count += 1
        if all_kmers_sum[k] == all_kmers[k]:
            print(f">{code[all_kmers[k]]}_{count}\n{k}")




if __name__ == '__main__':
    main()


