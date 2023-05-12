#!/usr/bin/env python3
import sys
from icecream import ic
if len(sys.argv) < 2:
    sys.exit(f"python3 {sys.argv[0]} all.hifi.rdna.variants.readbase.tab")



def addtodict3(thedict,key_a,key_b,key_c,val):
    # key_a = read id
    # key_b = ref id
    # key_c = ref' position
    # val = genotype
    if key_a in thedict:
        if key_b in thedict[key_a]:
            thedict[key_a][key_b].update({key_c:val})
        else:
            thedict[key_a].update({key_b:{key_c:val}})
    else:
        thedict.update({key_a:{key_b:{key_c:val}}})

def main():
    nucMatrix = {}
    ref_scafs = {}
    with open(sys.argv[1], 'r') as fh:
        for i in fh:
            # KY962518.1    640 G   T   m64144_211230_215102/7472316/ccs    G   SAME
            [scaf, pos, ref, alt, read, var, signal] =  i.strip().split("\t")
            genotype = '0'  # same
            if signal == 'DIFF':
                genotype = '2'
            if scaf not in ref_scafs:
                ref_scafs[scaf] = [pos, ]
            else:
                if pos not in ref_scafs[scaf]:
                    ref_scafs[scaf].append(pos)
            addtodict3(nucMatrix, read, scaf, pos, genotype)
    # generate output file header
    header = ['Reads', ]
    #header = []
    for s in sorted(ref_scafs):
        for p in sorted(ref_scafs[s]):
            header.append(s + ":" + p)
    print("\t".join(header))
    for r in nucMatrix:
        output = [r, ]
        for scaf in sorted(ref_scafs):
            for pos in ref_scafs[scaf]:
                if pos in nucMatrix[r][scaf]:
                    output.append(nucMatrix[r][scaf][pos])
                else:
                    output.append('0.5') # unknown 
        print("\t".join(output))


if __name__ == '__main__':
    main()

