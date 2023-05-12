#!/usr/bin/env python3
import sys

if len(sys.argv) < 2:
    sys.exit(f"python3 {sys.argv[0]} syri.vcf")


def main():
    with open(sys.argv[1], 'r') as fh:
        for i in fh:
            if i.startswith("#"):continue
            tmp = i.strip().split()
            # CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
            chrA = tmp[0]
            chrA_start = int(tmp[1])
            varid = tmp[2]
            ref_seq = tmp[3]
            alt_seq = tmp[4]
            tmp1 = tmp[7].split(";")
            chrA_end = tmp1[0].replace('END=', '')
            chrA_end = int(chrA_end)
            chrB = tmp1[1].replace('ChrB=', '')
            chrB_start = tmp1[2].replace('StartB=', '')
            chrB_end = tmp1[3].replace('EndB=', '')
            vartype = varid
            for c in "0123456789":
                vartype = vartype.replace(c, '')
            if vartype in ['INVAL', 'INVDPAL', 'INVTRAL', 'TRANSAL', 'SYNAL', 'NOTAL', 'DUPAL']:
                continue
            else:
                if vartype == 'SNP':
                    varlen = 1
                elif vartype == 'DEL':
                    varlen = len(ref_seq) - len(alt_seq)
                elif vartype == 'INS':
                    varlen = len(alt_seq) - len(ref_seq)
                else:
                    varlen = int(chrA_end) - int(chrA_start) + 1
                # vcf 2 bed. pos - 1
                chrB_start, chrB_end = int(chrB_start), int(chrB_end)
                print(f"{chrA}\t{chrA_start-1}\t{chrA_end-1}\t{vartype}\t{varid}\t{varlen}\t{chrB}\t{chrB_start-1}\t{chrB_end-1}\t{ref_seq}\t{alt_seq}")


if __name__ == '__main__':
    main()

