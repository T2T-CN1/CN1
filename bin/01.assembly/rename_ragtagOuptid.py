#!/usr/bin/env python3
import sys
if len(sys.argv) < 5:
    print(f"python3 {sys.argv[0]} ragtag.patch.ctg.agp ragtag.patch.rename.agp ragtag.patch.agp chr")
    exit()

def ctg_id(File):
    ctgs = {}
    with open(File, 'r') as fh:
        for i in fh:
            if i.startswith("#"):
                continue
            tmp = i.strip().split()
            ctgs[tmp[5]] = tmp[0]

    return ctgs

def qry_id(File):
    qrys = {}
    with open(File, 'r') as fh:
        for i in fh:
            if i.startswith("#"):
                continue
            tmp = i.strip().split()
            qrys[tmp[0]] = tmp[5]

    return qrys

def main():
    output = []
    ctgs = ctg_id(sys.argv[1])
    qrys = qry_id(sys.argv[2])
    out = open("ragtag.oriID.agp", 'w')
    out1 = open("ragtag.mapID.txt", 'w')
    with open(sys.argv[3], 'r') as fh:
        for i in fh:
            if i.startswith("#"):
                print(i.strip(), file=out)
            else:
                tmp = i.strip().split()
                if tmp[0].startswith("scf"):
                    new_id = sys.argv[4] + "_rg" + tmp[0][-3:]
                    pair = f"{tmp[0]}\t{new_id}"
                    if pair not in output:
                        output.append(pair)
                    tmp[0] = new_id

                elif tmp[0].startswith("seq"):
                    new_id = ctgs[tmp[0]]
                    pair = f"{tmp[0]}\t{new_id}"
                    if pair not in output:
                        output.append(pair)
                    tmp[0] = new_id
                if tmp[5].startswith("q"):
                    tmp[5] = qrys[tmp[5]]
                else:
                    tmp[5] = ctgs[tmp[5]]

                print("\t".join(tmp), file=out)
    out.close()
    for p in output:
        print(p, file=out1)
    out1.close()

if __name__ == '__main__':
    main()
