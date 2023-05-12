#/usr/bin/env python3
import sys
import os
from icecream import ic

import pyfastx

if len(sys.argv) < 3:
    sys.exit(f"python3 {sys.argv[0]} *.cent.bed original.fasta")

def main():
    centromeres = {}
    with open(sys.argv[1], 'r') as fh:
        for i in fh:
            tmp = i.strip().split()
            start = int(tmp[1])
            end = int(tmp[2])
            if tmp[0] not in centromeres:
                centromeres[tmp[0]] = [(start, end),]
            else:
                centromeres[tmp[0]].append((start, end))

    fa = pyfastx.Fastx(sys.argv[2])
    for name,seq,comment in fa:
        physic_end = len(seq.strip())
        if name in centromeres:
            # this chr has centromere, and i will break it into pieces
            index = 1
            current_start = 0
            for bnd in sorted(centromeres[name]):
                ic(bnd)
                start, end = bnd
                if index == 1:
                    # this is first one
                    # some chrs are acrocentric
                    if start == current_start:
                        current_start = end
                        continue
                print(f">{name}_{current_start}-{start}\n{seq[current_start:start]}")
                current_start = end
                index += 1
            if current_start < physic_end:
                print(f">{name}_{current_start}-{physic_end}\n{seq[current_start:]}")


if __name__ == '__main__':
    main()
