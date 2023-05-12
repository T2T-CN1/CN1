#!/usr/bin/env python3

import sys
import argparse
from icecream import ic


def getLen(length_file):
    lens = {}
    with open(length_file, 'r') as fh:
        for i in fh:
            tmp = i.strip().split()
            lens[tmp[0]] = int(tmp[1])
    return lens


def get_coveraged_length(bed_list):
    # merge bed to generate non-overlap blocks
    mapped_length = 0
    sort_bed_list = sorted(bed_list)

    low_est = sort_bed_list[0][0]
    high_est = sort_bed_list[0][1]

    for index, block in enumerate(sort_bed_list):
        low, high = block
        if high_est >= low:
            if high_est < high:
                high_est = high
        else:
            mapped_length += (high_est - low_est + 1)
            low_est, high_est = sort_bed_list[index]
    mapped_length += (high_est - low_est + 1)
    return mapped_length


def Add_records(dict, k, v):
    if k in dict.keys():
        dict[k].append(v)
    else:
        dict[k] = [v,]

def judge_orientation(orientation_f, orientation_r):
    results = {}
    qf = set(orientation_f.keys())
    qr = set(orientation_r.keys())
    qrys = set.union(qf, qr)
    for i in qrys:
        if i in qf and i not in qr:
            results[i] = "+"
        elif i in qf and i in qr:
            if orientation_f[i] >= orientation_r[i]:
                results[i] = "+"
            else:
                results[i] = "-"
        elif i not in qf and i in qr:
            results[i] = "-"
        else:
            sys.exit("something wrong")

    return results

def sortbyCoords(kary_que, ref_coverage):
    middle_position = {}
    for k in kary_que:
        sorted_blocks = sorted(ref_coverage[k],key=lambda x:(x[0],x[1]))
        very_begining = sorted_blocks[0][0]
        very_ending = sorted_blocks[-1][-1]
        middle = round((very_begining + very_ending ) / 2)
        middle_position[k] = middle

    sort_by_coords = sorted(middle_position.keys(), key=lambda k:middle_position[k])
    return sort_by_coords

def main(args):
    synteny = {}
    query_coverage = {}
    ref_coverage = {}
    scaf_length = {}
    orientation_f = {}
    orientation_r = {}

    with open(args.paf, 'r') as fh:
        # input file is in paf format
        for b in fh:
            tmp = b.strip().split()
            strand = tmp[4]
            ref_start = int(tmp[7])
            ref_end = int(tmp[8])
            qry_start = int(tmp[2])
            qry_end = int(tmp[3])
            ref_alen  = abs(ref_end - ref_start)
            qry_aln   = abs(qry_end - qry_start)
            identity  = int(tmp[9]) / int(tmp[10])
            ref_length = int(tmp[6])
            qry_length = int(tmp[1])
            ref_id = tmp[5]
            qry_id = tmp[0]
            scaf_length[ref_id] = ref_length
            scaf_length[qry_id] = qry_length
            # orientation against ref
            if strand == '+':
                if qry_id not in orientation_f:
                    orientation_f[qry_id] = qry_aln
                else:
                    orientation_f[qry_id] += qry_aln
            else:
                if qry_id not in orientation_r:
                    orientation_r[qry_id] = qry_aln
                else:
                    orientation_r[qry_id] += qry_aln

            # make all positive strand
            if qry_start > qry_end:
                qry_end, qry_start = qry_start, qry_end
            # calculate query's coverage before filtering
            if qry_id not in query_coverage:
                query_coverage[qry_id] = [[qry_start, qry_end],]
            else:
                query_coverage[qry_id].append([qry_start, qry_end])

            # prepare to sort query by it's mapped position in ref
            if qry_id not in ref_coverage:
                ref_coverage[qry_id] = [[ref_start, ref_end],]
            else:
                ref_coverage[qry_id].append([ref_start, ref_end])

            if args.filter:
                # round 1, filtering by scaffold length and identity, and work out query coverage
                if ref_length < args.minL or qry_length < args.minL:
                    continue
                elif ref_alen < args.minB or qry_aln < args.minB:
                    continue
                elif identity < args.identity:
                    continue
            record = [ref_id, ref_start, ref_end, qry_id, qry_start, qry_end]
            Add_records(synteny, ref_id, record)

    if args.filter:
        query_remove = []
        # round 2, filtering by query coverage
        for q in query_coverage:
            coveraged_length = get_coveraged_length(query_coverage[q])
            if coveraged_length / scaf_length[q] < args.rate:
                query_remove.append(q)
        # then remove coords of removed querys
        for r in synteny.keys():
            fixed_records = synteny[r]
            new_records = []
            for coord in fixed_records:
                if coord[3] not in query_remove:
                    new_records.append(coord)
            synteny[r] = new_records


    kary_que = []
    for x in synteny:
        for i in synteny[x]:
            if i[3] not in kary_que:
                kary_que.append(i[3])

    # default order is sorted by coordinates in ref;
    new_que = sortbyCoords(kary_que, ref_coverage)
    orientation_decision = judge_orientation(orientation_f, orientation_r)
    for i in new_que:
        print(f"{i}\t{scaf_length[i]}\t{orientation_decision[i]}")


if __name__ == '__main__':
    des = """

    deal synteny block file and generate color profile for circos,
    imput file should be six-column table, querry_scaf, q_start,
    q_end, ref_scaf, r_start, r_end, also you need give a table
    containing scaffold's length.

    yangchentao at genomics.cn, BGI.
    """
    parser = argparse.ArgumentParser(description=des,
                                    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-paf", type=str, required=True, metavar='<File>',
                        help="alignment result (paf format) generated by unimap")
    parser.add_argument("-order", type=str, required=False, metavar='<File>',
                       help="guide order of karyotype to show in circos")
    parser.add_argument("-filter", action="store_true",
                       help="filter records or not?")
    parser.add_argument("-minL", type=int, required=False, metavar='<INT>',
                        help="min scaffold length in link record, default=1000",
                        default=1000)
    parser.add_argument("-minB", type=int, required=False, metavar='<INT>',
                        help="min length of synteny block, default=500",
                       default=500)
    parser.add_argument("-id", type=float, required=False, metavar='<FLOAT>', dest="identity",
                       help="min identiy of synteny block, default=95", default=0.90)
    parser.add_argument("-rate", type=float, required=False, metavar='<FLOAT>',
                        help="remove these query which coverage is shorter than scaf_len * rate, default=0.1",
                       default=0.1)
    args =  parser.parse_args()
    main(args)
