#!/usr/bin/env python3
import sys

if len(sys.argv) < 2:
    sys.exit(f"print {sys.argv[0]} *.gfa")

def main():
    with open(sys.argv[1], 'r') as fh:
        for i in fh:
            if i[0] == 'S':
                tmp = i.strip().split()
                tmp[2] = "*"
                print("\t".join(tmp))
            else:
                print(i.strip())


if __name__ == '__main__':
    main()
