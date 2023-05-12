#!/usr/bin/env python3
import sys

if len(sys.argv) < 2:
    sys.exit(f"python3 {sys.argv[0]} *.vcf")

def make_dict(string):
    thedict = {}
    arr = string.split(";")
    for item in arr:
        if '=' in item:
            a, b = item.split("=")
            thedict[a] = b

    return thedict


def main():
    # PRECISE;SVTYPE=DEL;SVLEN=-35;END=822463;SUPPORT=64;COVERAGE=60,61,63,64,62;STRAND=+-;AF=1.000;STDEV_LEN=0.000;STDEV_POS=0.000   GT:GQ:DR:DV     1/1:60:0:64
    with open(sys.argv[1], 'r') as fh:
        for i in fh:
            if i.startswith("#"):
                continue
            tmp = i.strip().split()
            # #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE
            ref, pos, svid, ref_seq, alt_seq, qual, filt, info, form, gt_info = tmp
            if "BND" in svid:
                continue
            pos = int(pos)
            arrt = make_dict(info)
            svtype = arrt['SVTYPE']
            svlen = arrt['SVLEN']
            svlen = abs(int(svlen))
            if svtype == "DEL" or svtype == "INV":
                start = pos - 1
                end = start + svlen
            elif svtype == "INS":
                start = pos - 1
                end = start + 1
            print(f"{ref}\t{start}\t{end}\t{svtype}\t{svlen}\t{form}\t{gt_info}")


if __name__ == '__main__':
    main()

