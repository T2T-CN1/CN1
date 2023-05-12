import argparse
import sys
import textwrap
from copy import copy

parser = argparse.ArgumentParser(prog='fastaKit',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
Description
    deal with pag issues
'''))
parser.add_argument('-paf',  metavar='<file>', help='input paf file')
parser.add_argument('-l',  dest="scaf_list", metavar='<list>', help='scaffold list that needs to be inverted, scaf1,scaf2,...')
parser.add_argument('-f',  dest="scaf_file", metavar='<file>', help='scaffold list file that needs to be inverted')

if len(sys.argv) == 1:
    parser.print_help()
    parser.exit()

args = parser.parse_args()

def main():
    transset = []
    if args.scaf_list:
        transset = args.scaf_list.strip().split(",")
    elif args.scaf_file:
        with open(args.scaf_file, "r") as fh:
            for i in fh:
                transset.append(i.strip())
    elif args.scaf_list and args.scaf_file:
        sys.exit(f"can not using both records")
    else:
        sys.exit(f"you must give a list of scaffolds or a file")

    with open(args.paf, "r") as paf:
        for line in paf:
            lining = line.rstrip().split()
            if lining[0] in transset:
                length = int(lining[1])
                start = int(lining[2])
                end = int(lining[3])
                newend = length - start
                newstart = length - end
                lining[0] = "%s_R"%(lining[0])
                lining[2] = str(newstart)
                lining[3] = str(newend)
                if lining[4] == "-":
                    lining[4] = "+"
                elif lining[4] == "+":
                    lining[4] == "-"
            outstr = " ".join(lining)
            print(outstr)

if __name__ == "__main__":
        main()
