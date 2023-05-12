#!/usr/bin/env python3
import sys
import os
import subprocess
from icecream import ic

if len(sys.argv) < 3:
    sys.exit(f"python3 {sys.argv[0]} bam bed")


def main():
    deal_scafs = set()
    with open(sys.argv[2], 'r') as fh:
        for i in fh:
            tmp = i.strip().split("\t")
            scaf, start, end = tmp[0], tmp[1], tmp[2]
            if scaf not in deal_scafs:
                deal_scafs.add(scaf)
                if os.path.exists(scaf) == False:
                    os.mkdir(scaf)
            outpre = scaf + "_" + start + "-" + end
            cmd = f'samtools view -h -q 10 -bS -o {scaf}/{outpre}.bam {sys.argv[1]} {scaf}:{start}-{end}'
            cmd1 = f'samtools index {scaf}/{outpre}.bam'
            cmd2 = f'samtools depth {scaf}/{outpre}.bam > {scaf}/{outpre}.depth'
            subprocess.call(cmd, shell=True)
            subprocess.call(cmd1, shell=True)
            subprocess.call(cmd2, shell=True)


if __name__ == '__main__':
    main()

