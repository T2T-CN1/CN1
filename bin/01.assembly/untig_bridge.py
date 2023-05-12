#!/opt/homebrew/opt/python@3.9/bin/python3.9
import os
import sys
import time
import gzip
import pyfastx
import argparse
import subprocess
from icecream import ic
t = time.time()


###############################################################################
#####------------------------- parameters --------------------------------#####


description = """

Description

    combine two assemblies by overlap

"""

parser = argparse.ArgumentParser(
    prog="HiFiam_Verkko_bridge",
    description=description,
    formatter_class=argparse.RawTextHelpFormatter,
)

parser.add_argument('FASTA', metavar='<FILE>',
                    type=str,
                    help="input fastas separted by ,")

parser.add_argument('-id', metavar='<FLOAT>',
                    type=float, required=False,
                    dest='overlap_identity',
                    default=0.98,
                    help="the minimum identity allowed in overlapping region, default=0.98")

parser.add_argument('-overlap', metavar='<INT>',
                    type=int, required=False,
                    default=10000, dest='min_overlap',
                    help="the minimum overlapped length, default=10000")

parser.add_argument('-outpre', metavar='<STR>',
                    type=str, required=False,
                    default="assembly",
                    help="output preix, default=assembly")


###############################################################################
#####---------------------- program execution start ----------------------#####

def print_time(str):
    print(str + " " + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
# ----------------------------------------------------------------

def check_and_open_outhandle(file):
    if os.path.exists(file):
        print("WARRNING: " + file + " exists! now overwriting ...")
    else:
        print("[INFO]: " + "open file " + file + "...")
    out = open(file, 'w')
    return out


def complementation(sequence):
    # make a sequence complement #
    # replace function of string is too low!
    sequence = sequence.upper() ## [a bug fixed], reported by Wu Ping 20181129
    transtable = str.maketrans('ATCG-', 'TAGC-')
    sequence = sequence.translate(transtable)
    return sequence

def comp_rev(sequence):
    # make a sequence complement and reversed #
    sequence = complementation(sequence)
    return sequence[::-1]

def match(str1, str2):
    # ----count matched bases of two sequences----#
    matched = 0
    for base in range(len(str1)):
        if str1[base] == str2[base]:
            matched += 1
    identity = matched / len(str1)
    return identity


def generate_consensus(nameA, seqA, nameB, seqB, step):
    path = None
    # by default, seqA is longer than seqB
    lenA = len(seqA)
    lenB = len(seqB)
    if step <= lenB:
        consensus = seqA[0:-step] + seqB
        path = f">{nameA}:0:{lenA - step}>{nameB}:0:{lenB}"
    elif step > lenB and step <= lenA:
        consensus = seqA
        path = f">{nameA}:0:{lenA}"
    else:
        consensus = seqB[0:step-lenA] + seqA
        path = f">{nameB}:0:{step-lenA}>{nameA}:0:{lenA}"

    return consensus, path

def MergeOverlapSegment(longerSeq, shorterSeq):
    # part 1, let the short sequence be positive
    ##-----------anchoring overlap site--------#
    lenA = len(longerSeq)
    lenB = len(shorterSeq)
    results = []
    if args.min_overlap > lenB:
        sys.exit("min overlap length must be less than short seuqnce's length")
    else:
        for step in range(args.min_overlap, lenA+lenB-args.min_overlap):
            """
                        step
                          |---la---|
            ------------------------
                          ------------------
                          |---lb---|
                step = 10
                lenA = 32
                lenB = 15
            """
            if step <= lenB:
                la = longerSeq[-step:]
                lb = shorterSeq[0:step]
            elif step > lenB and step <= lenA:
                """
                  step
                -------------------------
                  ----------------
                """
                la = longerSeq[-step:lenB-step]
                lb = shorterSeq
            else:
                """
            step
                 |----la----|
                 ------------------------
            -----------------
                 |----lb----|
                step = 29
                lenA = 24
                lenB = 17
                """
                la = longerSeq[0:lenB-step]
                lb = shorterSeq[step-(lenA+lenB):]
            tmp_identity = match(la, lb)
            overlap_length = len(la)
            #ic(step, la, lb, tmp_identity)
            if tmp_identity == 1:
                results.append((step, 1, overlap_length))
                # find best result, so exit loop #
                break

            elif tmp_identity >= args.overlap_identity:
                results.append((step, tmp_identity, overlap_length))

        # find best overlaping result in all potenial positions
        candidates = sorted(results,key=lambda x:(x[1],x[2]))

        return candidates

# ------------------------assembly process--------------------------#
def Assembly(args):
    print_time("[INFO]: Assembling start:")


    assembly_result = args.outpre + "_assembly.fasta"
    fh_out = check_and_open_outhandle(assembly_result)
    fh_log = check_and_open_outhandle(args.outpre + "_assembly.log")


    fh_log.write("## overlaping identity = "
                 + str(args.overlap_identity)
                 + "\n")

    # --------------main-----------------------#
    tmp = args.FASTA.strip().split(",")
    fa1 = pyfastx.Fastx(tmp[0])
    fa2 = pyfastx.Fastx(tmp[1])

    for name1,seq1,comment1 in fa1:
        for name2,seq2,comment2 in fa2:
            seq1 = seq1.upper()
            seq2 = seq2.upper()
            len1 = len(seq1)
            len2 = len(seq2)
            if len1 >= len2:
                seqA, seqB = seq1, seq2
                nameA, nameB = name1, name2
            else:
                seqA, seqB = seq2, seq1
                nameA, nameB = name2, name1
            potiential_assembly_positive = MergeOverlapSegment(seqA, seqB)

            if len(potiential_assembly_positive) > 0:
                step, oid, overlap_len = potiential_assembly_positive[0]
                this_oid = oid * 100
                makeup_consensus, path = generate_consensus(nameA, seqA, nameB, seqB, step)
                len_makeup_consensus = len(makeup_consensus)
                fh_out.write(f">{path};len={len_makeup_consensus};oid={this_oid}%\n{makeup_consensus}")

                # output all results into log
                fh_log.write(f">{nameA}:+ {nameB}:+\n")
                for index,assembly in enumerate(potiential_assembly_positive):
                    s, oid, l = assembly
                    asm,path = generate_consensus(nameA, seqA, nameB, seqB, s)
                    asmlen = len(asm)
                    fh_log.write(f"[{index}] {path}, step={s}, overlap_identity={oid}, overlap_length={l}, assembly_len={asmlen}\n")

                # if check result, and ok so write into output_checked file #
            else:
                fh_log.write(f"#{nameA}:+ {nameB}:+ has no result!\n")
                # NEXT STEP, check another strand
                seqB = comp_rev(seqB)
                nameB = nameB + ":-"
                potiential_assembly_negative = MergeOverlapSegment(seqA, seqB)
                if len(potiential_assembly_negative) > 0:
                    step, oid, overlap_len = potiential_assembly_negative[0]
                    makeup_consensus, path = generate_consensus(nameA, seqA, nameB, seqB, step)
                    len_makeup_consensus = len(makeup_consensus)
                    this_oid = oid * 100
                    fh_out.write(f">{path};len={len_makeup_consensus};oid={this_oid}%\n{makeup_consensus}")

                    # output all results into log
                    fh_log.write(f">{nameA}:+ {nameB}\n")
                    for index,assembly in enumerate(potiential_assembly_negative):
                        s, oid, l = assembly
                        asm, path = generate_consensus(nameA, seqA, nameB, seqB, s)
                        asmlen = len(asm)
                        fh_log.write(f"[{index}] {path}, step={s}, overlap_identity={oid}, overlap_length={l}, assembly_len={asmlen}\n")

                    # if check result, and ok so write into output_checked file #
                else:
                    fh_log.write(f"#{nameA}:+ {nameB} has no result!\n")
    fh_out.close()

    print_time("[INFO]: Assembling done:")

if __name__ == '__main__':
    args = parser.parse_args()
    Assembly(args)
