import sys
import re
if len(sys.argv) < 3:
    print("Usage: python3 {} {} {}".format(sys.argv[0], "out.snps", "prefix"))
    exit()

sh = open(sys.argv[2] + ".snp", 'w')
ih = open(sys.argv[2] + ".indel", 'w')
pattern = re.compile(r'^\d+')
with open(sys.argv[1], 'r') as fh:
    for i in fh:
        i = i.strip()
        if pattern.match(i) != None:
            tmp = i.split("\t")
            if tmp[1] == "." or tmp[2] == ".":
                print(i, file=ih)
            else:
                print(i, file=sh)

sh.close()
ih.close()

