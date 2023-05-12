import sys

if len(sys.argv) < 2:
    sys.exit(f"python3 {sys.argv[0]} *.bed(collapse)")


with open(sys.argv[1], 'r') as fh:
    for i in fh:
        tmp = i.strip().split()
        if tmp[-1] == ".":
            tmp[-1] = "0"
        else:
            vars = tmp[-1].split(",")
            var_count = len(vars)
            tmp[-1] = str(var_count)
        print("\t".join(tmp))
