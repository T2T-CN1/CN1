#!/usr/bin/env python3
import sys
import re
#from icecream import ic

help="""
About the description of definition of all kinds of variants,
see https://github.com/schneebergerlab/syri/issues/67
DUP is for larger structural rearrangements, whereas CPG/CPL and TDM corresponds
to "local" structural variations identified within syntenic or structurally
rearranged regions. So, it is possible, for example, to have a local TDM within
a syntenic or a CPG within a translocation. Consider it like nested genomic
variation identification, where first the genomes are divided into syntenic
and structural rearrangements, followed by identification of structural
variations (indels, TDM, HDR) within these syntenic and structural rearrangements.
"""
if len(sys.argv) < 2:
    print("Note:\n" + help)
    sys.exit("Usage: python3 {} {}".format(sys.argv[0], "syri.out"))

def main():
    with open(sys.argv[1], 'r') as fh:
        for i in fh:
            tmp = i.strip().split("\t")
            vartype = tmp[-2]
            note = tmp[-1]
            if vartype.endswith("AL"):
                continue
            ref_id = tmp[0]
            ref_start = int(tmp[1])
            ref_end = int(tmp[2])
            que_id = tmp[5]
            que_start = int(tmp[6])
            que_end = int(tmp[7])
            # make query's location as '+' strand
            if que_end < que_start:
                que_start, que_end = que_end, que_start
                tmp[6], tmp[7] = tmp[7], tmp[6]

            if vartype in ["SNP", "SYN"]:
                continue
            elif vartype == "HDR":
                ref_sv_len = abs(ref_end - ref_start) + 1
                que_sv_len = abs(que_end - que_start) + 1
                sv_len = ref_sv_len - que_sv_len
                if sv_len > 0 and sv_len >= 50:
                    tmp[-2] = 'HDR|deletion'
                    print("\t".join(tmp[:-1]))
                elif sv_len < 0 and sv_len <= -50:
                    tmp[-2] = 'HDR|insertion'
                    print("\t".join(tmp[:-1]))

            # for SNP, I will not use results generated form syri, instead,
            # I use SNP from nucmer/dnadiff
            elif vartype == "INS":
                sv_len = abs(que_end - que_start) + 1
                if sv_len >= 50:
                    tmp[3] = '-'
                    tmp[4] = '-'
                    print("\t".join(tmp[:-1]))
            elif vartype == "DEL":
                # only keep large indels (>50 bp), small indels will generated from nucmer/dnadiff
                sv_len = abs(ref_end - ref_start) + 1
                if sv_len >= 50:
                    tmp[3] = '-'
                    tmp[4] = '-'
                    print("\t".join(tmp[:-1]))
            elif vartype in ['CPG', 'CPL']:
                # local copy number changes
                tmp[-2] = 'CNV|' + tmp[-2]
                print("\t".join(tmp[:-1]))
            elif vartype in ['DUP', 'INVDP'] and note in ['copyloss', 'copygain']:
                # copy number change of duplication, 
                tmp[-2] = 'CNV|' + note
                print("\t".join(tmp[:-1]))
            elif vartype == 'TDM':
                # Tandem contraction/expansion
                ref_var_len = abs(ref_end - ref_start) + 1
                que_var_len = abs(que_end - que_start) + 1
                # unbalanced tandem, so it's contraction or expansion
                if ref_var_len >= 50 or que_var_len >= 50:
                    if ref_var_len > que_var_len:
                        tmp[-2] = 'CNV|TDM_contraction'
                        print("\t".join(tmp[:-1]))
                    elif ref_var_len < que_var_len:
                        tmp[-2] = 'CNV|TDM_expansion'
                        print("\t".join(tmp[:-1]))
            else:
                # INV, INVTR, TRANS
                print("\t".join(tmp[:-1]))


if __name__ == '__main__':
    main()

