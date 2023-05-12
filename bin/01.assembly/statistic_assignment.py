#!/usr/bin/env python3
import sys

if len(sys.argv) < 2:
    sys.exit(f"print {sys.argv[0]} *.paf")


def cal_average(alns):
    alnrate = []
    identity = []
    for r, i in alns:
        alnrate.append(r)
        identity.append(i)

    average_alnrate = sum(alnrate) / len(alnrate)
    average_identity = sum(identity) / len(identity)
    return average_alnrate, average_identity



def main():
    mat_assign = {}
    pat_assign = {}
    scafs = set()
    mat_best_hit = {}
    pat_best_hit = {}
    with open(sys.argv[1], 'r') as fh:
        for i in fh:
            tmp = i.strip().split()
            scafs.add(tmp[0])
            hap = tmp[5].split("_")[0]
            align_rate = int(tmp[9]) / int(tmp[1])
            identity = int(tmp[9]) / int(tmp[10])

            if hap == 'mat':
                if tmp[0] not in mat_assign:
                    mat_assign[tmp[0]] = [(align_rate, identity),]
                else:
                    mat_assign[tmp[0]].append((align_rate, identity))
                if tmp[0] not in mat_best_hit:
                    mat_best_hit[tmp[0]] = identity
                else:
                    if identity > mat_best_hit[tmp[0]]:
                        mat_best_hit[tmp[0]] = identity
            elif hap == 'pat':
                if tmp[0] not in pat_assign:
                    pat_assign[tmp[0]] = [(align_rate, identity),]
                else:
                    pat_assign[tmp[0]].append((align_rate, identity))
                if tmp[0] not in pat_best_hit:
                    pat_best_hit[tmp[0]] = identity
                else:
                    if identity > pat_best_hit[tmp[0]]:
                        pat_best_hit[tmp[0]] = identity


    for scaf in scafs:
        if scaf in mat_assign:
            mat_mapped_average_alignrate, mat_mapped_average_identity = cal_average(mat_assign[scaf])
        else:
            mat_mapped_average_alignrate, mat_mapped_average_identity = 0, 0

        if scaf in pat_assign:
            pat_mapped_average_alignrate, pat_mapped_average_identity = cal_average(pat_assign[scaf])
        else:
            pat_mapped_average_alignrate, pat_mapped_average_identity = 0, 0
        if scaf in mat_best_hit:
            mat_best = mat_best_hit[scaf]
        else:
            mat_best = 'NA'
        if scaf in pat_best_hit:
            pat_best = pat_best_hit[scaf]
        else:
            pat_best = 'NA'
        print(f"{scaf}\tmat_best:{mat_best}\tmat_ave:{mat_mapped_average_alignrate}:{mat_mapped_average_identity}\tpat_best:{pat_best}\tpat_ave:{pat_mapped_average_alignrate}:{pat_mapped_average_identity}")

if __name__ == '__main__':
    main()

