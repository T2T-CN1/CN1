import sys
import re
if len(sys.argv) < 2:
    print("Usage: python3 {} <out.indel> ".format(sys.argv[0]))
    exit()


out = open(sys.argv[1] + ".report.txt", 'w')
print("#Chrom\tR_pos\tQ_pos\ttype\tsubtype\tlength\tR_base\tQ_base", file=out)
def read_block(fh):
    rc = ''
    qc = ''
    block = []
    for i in fh:
        tmp = i.strip().split("\t")
        ref_id = tmp[-2]
        query_id = tmp[-1]
        if rc == '':
            rc = ref_id
            block.append(i.strip())
        if qc == '':
            qc = query_id
        
        if ref_id != rc or query_id != qc:
            yield rc, qc, block
            block = []
            rc = ref_id
            qc = query_id
        else:
            block.append(i.strip())

    yield rc, qc, block



with open(sys.argv[1], 'r') as fh:
    for t in read_block(fh):
        var_insert = {}
        var_delete = {}
        pos = []
        kind = {}
        pos_match = {}
        rc, qc, block = t
        for b in block:
            tmp = b.split("\t")
            if tmp[1] == ".":
                kind[tmp[0]] = "i"
                # insertion compared to Ref
                if tmp[0] not in var_insert.keys():
                    pos.append(tmp[0])
                    var_insert[tmp[0]] = tmp[2]
                    pos_match[tmp[0]] = tmp[3]
                else:
                    var_insert[tmp[0]] += tmp[2]
            elif tmp[2] == ".":
                kind[tmp[0]] = "d"
                # deletion compared to Ref
                if tmp[3] not in var_delete.keys():
                    pos.append(tmp[0])
                    var_delete[tmp[3]] = tmp[1]
                    pos_match[tmp[0]] = tmp[3]
                else:
                    var_delete[tmp[3]] += tmp[1]
            else:
                print("this is not a indel")

        for p in pos:
            k = kind[p]
            if k == "i":
                r_pos = p
                q_pos = pos_match[p]
                indel_len = len(var_insert[r_pos])
                if indel_len == 1:
                    t = "." + var_insert[r_pos]
                    print("{}\t{}\t{}_\t{}\t{}\t{}\t{}\t{}\t{}".format(rc, qc, r_pos, q_pos, "ins", t,
                                                  1, ".", var_insert[r_pos]), file=out)
                else:
                    q_pos1 = int(q_pos) + indel_len - 1
                    t = ".S"
                    print("{}\t{}\t{}_\t{}\t{}\t{}\t{}\t{}\t{}".format(rc, qc, r_pos, q_pos + "-" + str(q_pos1),
                                                          "ins", t, indel_len, ".",
                                                          var_insert[r_pos]), file=out)
            else:
                r_pos = p
                q_pos = pos_match[p]
                indel_len = len(var_delete[q_pos])
                if indel_len == 1:
                    t = var_delete[q_pos] + "."
                    print("{}\t{}\t{}\t{}_\t{}\t{}\t{}\t{}\t{}".format(rc, qc, r_pos, q_pos, "del", t,
                                                  1, var_delete[q_pos], "."), file=out)
                else:
                    r_pos1 = int(r_pos) + indel_len - 1
                    t = "S."
                    print("{}\t{}\t{}\t{}_\t{}\t{}\t{}\t{}\t{}".format(rc, qc, r_pos + "-" + str(r_pos1),
                                                              q_pos, "del", t, indel_len,
                                                              var_delete[q_pos], "."), file=out)

out.close()
