#!/usr/bin/env python3
import sys

if len(sys.argv) < 2:
    sys.exit(f"python3 {sys.argv[0]} sv.vcf")

def makeup_into(raw, end):
    tmp = raw.split(";")
    a = ";".join(tmp[:3])
    b = "END=" + str(end)
    c = ";".join(tmp[4:])
    return a + ";" + b + ";" + c

def main():
    with open(sys.argv[1], 'r') as fh:
        for i in fh:
            if i.startswith("#"):
                print(i.strip())
            else:
                tmp = i.strip().split()
                if 'DEL' in tmp[2] or 'INS' in tmp[2]:
                    #this is a deletion, which we need to revise
                    ref_seq = tmp[3]
                    start = int(tmp[1])
                    ref_seq_len = len(ref_seq)
                    real_end = ref_seq_len - 1 + start
                    tmp[7] = makeup_into(tmp[7], real_end)
                    print("\t".join(tmp))
                else:
                    print(i.strip())

if __name__ == '__main__':
    main()


