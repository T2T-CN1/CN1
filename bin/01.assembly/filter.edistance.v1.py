#!/usr/bin/python
import pysam
import re
from sys import argv


path_in = argv[1]
path_out = argv[2]
samfile = pysam.AlignmentFile(path_in, "rb")
fo=open(path_out, 'w')
for line in samfile:
    count=0
    mapping=line.cigarstring
    ref=line.reference_name
    if (ref != None):
        edistance=line.get_tag('NM')
        if mapping == "201M" and edistance <=1:
                fo.write('\t'.join([line.qname,line.reference_name,str(edistance),mapping])+"\n")
        if edistance == 0:
            match=re.findall(r'(\d+M)',mapping)
            for i in range(0, len(match)):
                if(match[i] == "200M"):
                    fo.write('\t'.join([line.qname,line.reference_name,str(edistance),mapping])+"\n")
        if edistance == 1:
            match=re.findall(r'(\d+M)',mapping)
            for i in range(0, len(match)):
                match[i]=match[i].split('M')[0]
                count+=int(match[i])
            if count == 200:
                fo.write('\t'.join([line.qname,line.reference_name,str(edistance),mapping])+"\n")
samfile.close()
fo.close()
