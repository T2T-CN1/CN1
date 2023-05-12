#!/usr/bin/env python3

import sys
#from icecream import ic


def getLen(length_file):
    lens = {}
    with open(length_file, 'r') as fh:
        for i in fh:
            tmp = i.strip().split()
            lens[tmp[0]] = int(tmp[1])
    return lens


def get_coveraged_length(alns):
    # extract bed list form alns
    bed_list = alns

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
            mapped_length += (high_est - low_est )
            low_est, high_est = sort_bed_list[index]
    mapped_length += (high_est - low_est )
    return mapped_length


def Add_records(dict, ka, kb, val):
    if ka in dict:
        if kb in dict[ka]:
            dict[ka][kb].append(val)
        else:
            dict[ka][kb] = [val,]
    else:
        dict.update({ka: {kb: [val,]}})

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

def getaverage_identity(alns):
    # [ref_id, ref_start, ref_end, qry_id, qry_start, qry_end, identity]
    tmp = []
    for a in alns:
        tmp.append(float(a[-1]))
    average = "{:.5f}".format(sum(tmp) / len(alns))
    return average

def update_ori(mapping_ori, ref_id, qry_id, qry_aln, strand):
    if ref_id not in mapping_ori:
        mapping_ori.update({ref_id:{qry_id:(qry_aln, strand)}})
    else:
        if qry_id not in mapping_ori[ref_id]:
            mapping_ori[ref_id].update({qry_id:(qry_aln, strand)})
        else:
            last_record_len, last_record_ori = mapping_ori[ref_id][qry_id]
            if qry_aln > last_record_len:
                mapping_ori[ref_id][qry_id] = (qry_aln, strand)


def main():
    synteny = {}
    query_coverage = {}
    ref_coverage = {}
    scaf_length = {}
    mapping_ori = {}

    with open(sys.argv[1], 'r') as fh:
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

            if ref_id in synteny:
                synteny[ref_id].append([ref_start, ref_end])
            else:
                synteny[ref_id] = [[ref_start, ref_end], ]

    all_aligned = 0
    all_reflen = 0
    for r in synteny:
        mapping_results = {}
        alns = synteny[r]
        nonoverlap_ref_aligned = get_coveraged_length(alns)
        all_aligned += nonoverlap_ref_aligned
        all_reflen += scaf_length[r]
        alignrate = "{:.5f}".format(nonoverlap_ref_aligned / scaf_length[r])
        print(f"{r}\t{scaf_length[r]}\t{nonoverlap_ref_aligned}\t{alignrate}")
    overall_rate = "{:.5f}".format(all_aligned / all_reflen)
    print(f"overall\t{all_reflen}\t{all_aligned}\t{overall_rate}")


if __name__ == '__main__':
    des = """

    deal synteny block file and generate color profile for circos,
    imput file should be six-column table, querry_scaf, q_start,
    q_end, ref_scaf, r_start, r_end, also you need give a table
    containing scaffold's length.

    yangchentao at genomics.cn, BGI.
    """
    if len(sys.argv) < 2:
        sys.exit(f"python3 {sys.argv[0]} *.paf")
    else:
        main()
