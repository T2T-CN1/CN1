#!/usr/bin/env python3

import sys
import pysam
from collections import defaultdict

if len(sys.argv) != 3:
	sys.exit('python3 %s <target.bed> <bam>' % (sys.argv[0]))

bedFile = sys.argv[1]
bamFile = sys.argv[2]
samfile = pysam.AlignmentFile(bamFile, 'rb')

with open(bedFile) as f:
	for line in f:
		line = line.rstrip()
		tmp = line.split('\t')
		scaf = tmp[0]
		bg = int(tmp[1])
		ed = int(tmp[2])

		tmpDict = defaultdict(list)
		seqDict = {}
		for pileupcolumn in samfile.pileup(scaf, bg, bg+1):
			if pileupcolumn.pos == bg:    #only look at the speicific site
				for pileupread in pileupcolumn.pileups:
					if pileupread.indel == 0 and not pileupread.is_del and not pileupread.is_refskip:
						qname = pileupread.alignment.query_name
						qbg = pileupread.query_position
						tmpDict[qname].append(qbg)
						seqDict[qname] = pileupread.alignment.query_sequence

		for pileupcolumn in samfile.pileup(scaf, ed-1, ed):
			if pileupcolumn.pos == ed-1:    #only look at the speicific site
				for pileupread in pileupcolumn.pileups:
					if pileupread.indel == 0 and not pileupread.is_del and not pileupread.is_refskip:
						qname = pileupread.alignment.query_name
						qed = pileupread.query_position
						tmpDict[qname].append(qed)
						seqDict[qname] = pileupread.alignment.query_sequence
		
		for qname in tmpDict:
			l = len(tmpDict[qname])
			if l == 2:
				#pos = tmpDict[qname][0]
				seq = seqDict[qname][tmpDict[qname][0]:tmpDict[qname][1]+1]
				print('%s\t%i\t%i\t%s\t%i\t%i\t%s' % (scaf, bg, ed, qname, tmpDict[qname][0], tmpDict[qname][1]+1, seq))
			
			
