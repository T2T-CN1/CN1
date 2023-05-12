import pyfastx
import sys
if len(sys.argv) < 4:
    sys.exit(f"python3 {sys.argv[0]} *.fasta kmer_size step")

fa = pyfastx.Fastx(sys.argv[1])
kmer_size = int(sys.argv[2])
step = int(sys.argv[3])

uniq_kmers = {}
kmers = []

for name,seq,comment in fa:
    seqlen = len(seq)
    start = 0
    while start < seqlen - kmer_size:
        kmer = seq[start:start+kmer_size]
        real_start = start + 1
        real_end = real_start + kmer_size - 1q
        info = f"{name}:{real_start}-{real_end}"
        if kmer not in uniq_kmers:
            uniq_kmers[kmer] = [info,]
        else:
            uniq_kmers[kmer].append(info)
        start += step



for index,k in enumerate(sorted(uniq_kmers,key=lambda x:len(uniq_kmers[x]), reverse=True)):
    pos = ",".join(uniq_kmers[k])
    count = index + 1
    print(f"k_{count}\t{k}\t{pos}")
