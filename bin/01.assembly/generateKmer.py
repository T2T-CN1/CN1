import pyfastx
import sys
if len(sys.argv) < 3:
    sys.exit(f"python3 {sys.argv[0]} *.fasta kmer_size [uniq]?")

fa = pyfastx.Fastx(sys.argv[1])
start = 0
step = int(sys.argv[2])
uniq_kmers = []
kmers = []
for name,seq,comment in fa:
    seqlen = len(seq)
    while start < seqlen - step:
        kmer = seq[start:start+step]
        if kmer not in uniq_kmers:
            uniq_kmers.append(kmer)
        else:
            uniq_kmers.remove(kmer) # remove kmer that appears more than once
        kmers.append([start, kmer],)
        start += 1

if len(sys.argv) > 3 and sys.argv[3] == "uniq":
    for pos, kmer in kmers:
        if kmer in uniq_kmers:
            print(f"{pos+1}\t{kmer}")
else:
    for pos, kmer in kmers:
        print(f"{pos+1}\t{kmer}")
