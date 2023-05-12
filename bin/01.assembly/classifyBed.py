#!/usr/bin/env python3
import sys

def main():
    total = 0
    class_len = {}
    classification = {}
    with open(sys.argv[1], 'r') as fh:
        for i in fh:
            tmp = i.strip().split()
            if len(tmp) < 4:
                sys.exit("the bed file requires the fouth column")
            chrom, start, end, kind = tmp
            start = int(start)
            end = int(end)
            tmp_len = end - start
            if kind in class_len:
                class_len[kind] += tmp_len
                classification[kind].append(chrom+":"+str(start)+":"+str(end))
            else:
                class_len[kind] = tmp_len
                classification[kind] = [chrom+":"+str(start)+"-"+str(end),]

    for c in sorted(class_len.keys()):
        regions = ",".join(classification[c])
        print(f"{c}\t{class_len[c]}\t{regions}")

if __name__ == '__main__':
    main()
