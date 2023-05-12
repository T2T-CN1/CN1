#!/usr/bin/env python3
import sys
if len(sys.argv) < 3:
    sys.exit("python3 {} <scaf.fa> <*.agp>".format(sys.argv[0]))


def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, "".join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, "".join(seq))


def revcom_transtable(sequence):
    # make a sequence complement #
    transtable = str.maketrans('ATCGatcg', 'TAGCtagc')
    sequence = sequence.translate(transtable)
    return sequence[::-1]

def wrapseq(seq):
    l = len(seq)
    for i in range(0,l,80):
        if i+80 <= l:
            print(seq[i:i+80])
        else:
            print(seq[i:])

def main():
    scafs = {}
    need_curated_scafs = []
    with open(sys.argv[1], 'r') as fa:
        for name, seq in read_fasta(fa):
            name = name.replace(">", "").split()[0]
            scafs[name] = seq
    # need curated scaffolds AGP
    with open(sys.argv[2], 'r') as fh:
        data = {}
        superscafs = set()
        for i in fh:
            if i.startswith("#"):continue
            if len(i) == 0:continue
            tmp = i.strip().split()
            if len(tmp) != 9:
                sys.exit("bad AGP format of {}".format(sys.argv[2]))
            superscafs.add(tmp[0])
            if tmp[0] in data.keys():
                data[tmp[0]].append(tmp)
            else:
                data[tmp[0]] = [tmp,]
    # need curated scaffolds
    for s in list(superscafs):
        need_curated_scafs.append(s.split(".")[0]) # the scaffold id is like: Super-Scaffold_14.cur1, "Super-Scaffold_14" is original id
        makeup_superscaf = ""
        for r in data[s]:   # each subseq in super scaf
            pos_s = int(r[1])
            pos_e = int(r[2])
            region_size = pos_e - pos_s + 1
            comp_type = r[4]
            if comp_type != 'N':
                subscaf = r[5]
                subscaf_beg = int(r[6])
                subscaf_end = int(r[7])
                filled_size = subscaf_end - subscaf_beg + 1
                strand = r[8]
                if region_size != filled_size:
                    sys.exit("wrong position of {}".format(sys.argv[2]))
            else:
                gap_size = int(r[5])
                strand = '+' # meanless
                if region_size != gap_size:
                    sys.exit("wrong position of {}".format(sys.argv[2]))

            if strand == "+":
                if comp_type == 'W':
                    makeup_superscaf += scafs[subscaf][subscaf_beg-1:subscaf_end]
                elif comp_type == 'N':
                    makeup_superscaf += 'N'*gap_size
                else:
                    sys.exit("this type {} is not accepted".format(comp_type))
            else:
                if comp_type == 'W':
                    makeup_superscaf += revcom_transtable(scafs[subscaf][subscaf_beg-1:subscaf_end])
                elif comp_type == 'N':
                    makeup_superscaf += 'N'*gap_size
                else:
                    sys.exit("this type {} is not accepted".format(comp_type))
        print(">" + s)
        wrapseq(makeup_superscaf)

    # don't need to curate
    """
    for name in scafs:
        if name not in need_curated_scafs:
            print(">" + name)
            wrapseq(scafs[name])
    """

if __name__ == '__main__':
    main()
