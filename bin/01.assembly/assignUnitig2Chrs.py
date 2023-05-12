#!/usr/bin/env python3
import sys
import os
import pyfastx

if len(sys.argv) < 3:
    print("Usage: python3 {} <asm2ref.paf.stat> <assembly.fasta>\n".format(sys.argv[0]))
    exit()
#elif sys.argv[2] != "paf" and sys.argv[2] != "fmash":
#    print("only recognize paf or fmash (mashmap like) file format")


def reverseComplement(sequence):
    sequence = sequence.upper()
    transtable = str.maketrans('ATCGatcg', 'TAGCtagc')
    sequence = sequence.translate(transtable)
    return sequence[::-1]

fasta = {}
fa = pyfastx.Fastx(sys.argv[2])
for name,seq,comment in fa:
    fasta[name] = seq

records = {}

with open(sys.argv[1], 'r') as fh:
    for i in fh:
        tmp = i.strip().split('\t')
        unitig = tmp[0]
        best_hit = tmp[2]
        ref = best_hit.split()[0]
        if ref not in records:
            records[ref] = [[unitig, best_hit],]
        else:
            records[ref].append([unitig, best_hit])

if not os.path.exists("splitbyChrs"):
    os.mkdir("splitbyChrs")


for k in sorted(records.keys()):
    outdir = "splitbyChrs/" + k
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    aln = open(outdir + "/" + k + ".besthit" , 'w')
    fh = open(outdir + "/" + k + ".forward.fa", 'w')
    for i in records[k]:
        aln.write("\t".join(i))
        aln.write('\n')
        unitig, besthit = i
        strand = besthit.split()[-1]
        if strand == '-':
            seq = reverseComplement(fasta[unitig])
            unitig = unitig + '_R'
        else:
            seq = fasta[unitig]
        fh.write(f">{unitig}\n{seq}\n")
    aln.close()
    fh.close()
