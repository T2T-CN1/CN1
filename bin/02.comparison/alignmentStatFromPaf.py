#!/usr/bin/env python3

import sys
#from icecream import ic

import gzip
import os
def smart_open(file, opera):
    if opera == 'r':
        if os.path.exists(file) ==False:
            print("Can not open file {}".format(file))
            exit()
        else:
            if file.endswith(".gz"):
                out = gzip.open(file, 'rt')
            else:
                out = open(file, 'r')
    elif opera == 'w':
        if file.endswith(".gz"):
            out = gzip.open(file, 'wt')
        else:
            out = open(file, 'w')
    return out

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


def parser_blast(file):
    pre = ""
    target_match = []
    with smart_open(file, 'r') as fh:
        while True:
            line = fh.readline().strip()
            if not line:
                break
            tmp = line.strip().split("\t")
            if pre == "":
                target_match = [line,]
                pre = tmp[0]
            elif tmp[0] == pre:
                target_match.append(line)
            else:
                yield target_match
                target_match = [line,]
                pre = tmp[0]
        yield target_match

def main():
    input_paf = sys.argv[1]

    # input file is in paf format
    for records in parser_blast(input_paf):
        synteny = {}
        query_coverage = {}
        ref_coverage = {}
        scaf_length = {}
        mapping_ori = {}
        for b in records:
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

            # record mapping ori
            update_ori(mapping_ori, ref_id, qry_id, qry_aln, strand)
            # make all positive strand
            if qry_start > qry_end:
                qry_end, qry_start = qry_start, qry_end

            record = [ref_id, ref_start, ref_end, qry_id, qry_start, qry_end, identity]
            Add_records(synteny, qry_id, ref_id, record)


        for q in synteny:
            mapping_results = {}
            for r in synteny[q]:
                alns = synteny[q][r]
                average_identity = getaverage_identity(alns)
                qry_bed_list = [ [a[4],a[5]] for a in alns ]
                ref_bed_list = [ [a[1],a[2]] for a in alns ]
                nonoverlap_qry_aligned = get_coveraged_length(qry_bed_list)
                nonoverlap_ref_aligned = get_coveraged_length(ref_bed_list)
                qry_alignrate = "{:.5f}".format(nonoverlap_qry_aligned / scaf_length[q])
                ref_alignrate = "{:.5f}".format(nonoverlap_ref_aligned / scaf_length[r])
                # filtering low coverage results
                if float(qry_alignrate) < 0.4 or float(average_identity) < 0.2:
                    continue
                ori = mapping_ori[r][q][-1]
                info = f"{r} {scaf_length[r]} {nonoverlap_ref_aligned} {ref_alignrate} | {q} {scaf_length[q]} {nonoverlap_qry_aligned} {qry_alignrate} | {average_identity} {ori}"
                mapping_results[info] =  float(average_identity)
            # sorting mapping results by score (identity * alignrate)
            output = sorted(mapping_results, key=lambda k: (mapping_results[k], k), reverse=True)
            #stat = "\t".join(output)
            #print(f"{q}\t{scaf_length[q]}\t{stat}")
            for match in output:
                print(f"{match}")


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

