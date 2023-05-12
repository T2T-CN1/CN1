#!/usr/bin/env python3

import sys
import re
#import pyfastx
#from icecream import ic

if len(sys.argv) < 3:
    sys.exit(f"python3 {sys.argv[0]} merge_nodups.txt output.bed")

def sum_strs(strings):
    total = 0
    for s in strings:
        total += int(s)
    return total

def main():
    print("converting the format to bed and checking invalid symbols(:)...")
    output = sys.argv[2]
    bed_out = open(output, 'w')
    with open(sys.argv[1], 'r') as fh:
        for i in fh:
            tmp = i.split()
            str1, chr1, pos1, frag1, str2, chr2, pos2, frag2, mapq1, cigar1, sequence1, mapq2, cigar2, sequence2, readname1, readname2 = tmp
            if chr1 == chr2:
                matches1 = re.findall('(\d+)[MID]', cigar1)
                matches2 = re.findall('(\d+)[MID]', cigar2)
                #ic(matches1, matches2)
                bed_pos1 = int(pos1) - 1
                bed_pos2 = int(pos2) - 1
                end1 = bed_pos1 + sum_strs(matches1)
                end2 = bed_pos2 + sum_strs(matches2)
                if str1 == '0':
                    strand1 = '+'
                else:
                    strand1 = '-'

                if str2 == '0':
                    strand2 = '+'
                else:
                    strand2 = '-'
                print(f"{chr1}\t{bed_pos1}\t{end1}\t{readname1}/1\t{mapq1}\t{strand1}", file=bed_out)
                print(f"{chr2}\t{bed_pos2}\t{end2}\t{readname2}/2\t{mapq2}\t{strand2}", file=bed_out)

        bed_out.close()

"""
        fa = pyfastx.Fastx(sys.argv[2])
        fasta_out = open(sys.argv[2] + ".format_checked.fasta", 'w')
        for name,seq,comment in fa:
            if ":" in name:
                name = name.replace(":", "_")
                print(f">{name}\n{seq}", file=fasta_out)
        fasta_out.close()
"""


if __name__ == '__main__':
    main()

