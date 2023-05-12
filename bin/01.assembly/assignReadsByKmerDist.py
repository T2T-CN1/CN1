#/usr/bin/env python3

import sys
from icecream import ic

if len(sys.argv) < 2:
    sys.exit(f"python3 {sys.argv[0]} *.filter.edit1.txt")

def add2Ddict(thedict, ka, kb, val):
    if ka in thedict:
        if kb in thedict[ka]:
            thedict[ka].update({kb: val+thedict[ka][kb]})
        else:
            thedict[ka].update({kb: val})
    else:
        thedict.update({ka: {kb: val}})



def main():
    perfect_match_score = 5 # 201M
    notbad_match_score = 4 # 1D200M...
    kmer_matrix = {}
    with open(sys.argv[1], 'r') as fh:
        for i in fh:
            tmp = i.strip().split()
            kmer, readid, edit_distance, cigar = tmp
            chrom = kmer.split("_")[0]
            if int(edit_distance) == 0:
                add2Ddict(kmer_matrix, readid, chrom, perfect_match_score)
            elif int(edit_distance) == 1:
                add2Ddict(kmer_matrix, readid, chrom, notbad_match_score)

    for r in kmer_matrix:

        sorted_chrs = sorted(kmer_matrix[r], key=lambda a: (kmer_matrix[r][a], a), reverse=True)
        kmer_dist = []
        for k in sorted_chrs:
            kmer_dist.append(f"{k}:{kmer_matrix[r][k]}")
        tmp = "\t".join(kmer_dist)
        print(f"{r}\t{tmp}")


if __name__ == '__main__':
    main()




