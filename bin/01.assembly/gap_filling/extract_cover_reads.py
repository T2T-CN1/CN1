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

		seqDict = {}
		
		for pileupcolumn in samfile.pileup(scaf, bg, ed):
			if pileupcolumn.pos >= bg and pileupcolumn.pos < ed:
				for pileupread in pileupcolumn.pileups:
					if pileupread.indel == 0 and not pileupread.is_del and not pileupread.is_refskip:
						qname = pileupread.alignment.query_name	
						qpos = pileupread.query_position
						if qname not in seqDict:
							seqDict[qname] = [qpos, qpos+1]
						else:
							if qpos < seqDict[qname][0]: seqDict[qname][0] = qpos
							if qpos+1 > seqDict[qname][1]: seqDict[qname][1] = qpos+1
			
		for read in samfile.fetch(scaf, bg, ed):
			if read.query_name in seqDict:
				start, end = seqDict[read.query_name][0:2]
				print('%s\t%i\t%i\t%s\t%i\t%i\t%s' % (scaf, bg, ed, read.query_name, start, end, read.query_sequence[start:end]))
				
