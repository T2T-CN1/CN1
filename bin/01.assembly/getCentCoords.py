#!/usr/bin/env python3
import sys

if len(sys.argv) < 3:
    sys.exit(f"print {sys.argv[0]} *.paf cent.bed.txt")


def shift_pos(target, start, qry_start, qry_end, strand):
    if strand == "+":
        shift = target - start
        return qry_start + shift
    else:
        shift = target - start
        return qry_end - shift


def main():
    with open(sys.argv[2], 'r') as ch:
        tmp = ch.readline().strip().split()
        start = int(tmp[1])
        end = int(tmp[2])
    id1 = ""
    id2 = ""
    find1 = 0
    find2 = 0
    with open(sys.argv[1], 'r') as fh:
        for i in fh:
            tmp = i.strip().split()
            ref_start, ref_end = int(tmp[7]), int(tmp[8])
            qry_start, qry_end = int(tmp[2]), int(tmp[3])
            strand = tmp[4]
            if start > ref_start and start < ref_end:
                find1 = 1
                qry_cent_start = shift_pos(start, ref_start, qry_start, qry_end, strand)
                id1 = tmp[0]
            if end > ref_start and end < ref_end:
                find2 = 1
                qry_cent_end = shift_pos(end, ref_start, qry_start, qry_end, strand)
                id2 = tmp[0]
    if find1 and find2:
        if id1 == id2:
            print(f"{id1}\t{qry_cent_start}\t{qry_cent_end}")


if __name__ == '__main__':
    main()

