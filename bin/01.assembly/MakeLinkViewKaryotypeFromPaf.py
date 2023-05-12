#!/usr/bin/env python3
import sys
if len(sys.argv) < 2:
    sys.exit(f"python3 {sys.argv[0]} *.paf")

lens = {}
order = []
with open(sys.argv[1], 'r') as fh:
    for i in fh:
        [qry, qry_len, ref, ref_len] = i.strip().split()
        if qry not in lens:
            lens[qry] = qry_len
        if qry not in order:
            order.append(qry)

top  = []
for q in order:
    top.append(f"{q}:0:{lens[q]}")
print(" ".join(top))
print(f"{ref}:0:{ref_len}")
