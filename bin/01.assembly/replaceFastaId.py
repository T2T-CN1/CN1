import sys
def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name:
        yield (name, ''.join(seq))


if len(sys.argv) < 3:
    print("python3 {} {} fasta".format(sys.argv[0], "mapID"))
    exit()

records = {}
with open(sys.argv[1], 'r') as fh:
    for i in fh:
        tmp = i.strip().split()
        records[tmp[0]] = tmp[1]
with open(sys.argv[2], 'r') as fh:
    for i in read_fasta(fh):
        name, seq = i
        name = name.replace(">", "")
        if name in records:
            print(">" + records[name] + "\n" + seq)
        else:
            print(">" + name + "\n" + seq)


