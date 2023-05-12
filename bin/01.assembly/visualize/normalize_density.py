import sys
import math

with open(sys.argv[1], 'r') as fh:
    for i in fh:
        tmp = i.strip().split()
        freq = int(tmp[3])
        if freq == 0:
            print(i.strip())
        else:
            tmp[3] = str(math.log10(freq))
            print("\t".join(tmp))
