#!/usr/bin/env python3

import sys
from icecream import ic
import os

if len(sys.argv) < 2:
    sys.exit(f"python3 {sys.argv[0]} File1 File2 ..[File3]..")

def extract_name(filename):
    name = os.path.basename(filename)
    sample = name.replace(".bed", "")
    return sample


def main():
    # read the first file get each path in pangenome
    matrix = []
    header = ["#chrom", "start", "end", "pathIN", "pathOUT", ]
    with open(sys.argv[1], 'r') as fh:
        sample = extract_name(sys.argv[1])
        header.append(sample)
        for i in fh:
            tmp = i.rstrip().split()
            chrom, start, end, pathIN, pathOUT, pathINFO = tmp
            matrix.append(tmp,)


    # remain files
    appending_files = sys.argv[2:]
    for f in appending_files:
        with open(f, 'r') as fh:
            sample = extract_name(f)
            header.append(sample)
            index = 0
            for i in fh:
                tmp = i.rstrip().split()
                chrom, start, end, pathIN, pathOUT, pathINFO = tmp
                #ic(pathINFO)
                matrix[index].append(pathINFO)
                index  += 1


    print("\t".join(header))
    for item in matrix:
        print("\t".join(item))


if __name__ == '__main__':
    main()


