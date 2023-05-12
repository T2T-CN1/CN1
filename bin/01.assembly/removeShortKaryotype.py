#!/usr/bin/env python3
import sys

if len(sys.argv) < 3:
    sys.exit(f"python3 {sys.argv[0]} karyotype.txt outpre")


def main():
    outpre = sys.argv[2]
    out = open(outpre + ".short.karyotype.txt", 'w')
    long = open(outpre + ".1Mb.karyotype.txt", 'w')
    with open(sys.argv[1], 'r') as fh:
        for i in fh:
            new_list = []
            tmp = i.strip().split()
            for k in tmp:
                tmp1 = k.split(":")
                contig_len = int(tmp1[2])
                if contig_len > 1000000:
                    new_list.append(k)
                else:
                    out.write(tmp1[0] + "\t" + tmp1[2] + "\n")
            long_out = "\t".join(new_list)
            print(long_out, file=long)

    out.close()
if __name__ == '__main__':
    main()

