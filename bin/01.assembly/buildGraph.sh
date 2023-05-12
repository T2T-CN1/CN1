#!/bin/sh
set -e

#  Figure out what reads we're going to use.  If correction is enabled, the
#  reads are in a single file that is a relative path away from us; but if it
#  is not enabled, the reads are in possibly multiple files that are at an
#  absolute path.
#

if [[ $# != 3  ]] ; then
	echo "Usage : $0 hifi.fa kmer_size outpre"
		exit 1
fi

iopt=$1
kmer=$2
outpre=$3
cpu=8


echo "Building graph with reads:"
echo 

#  Build the graph.
#    (rule build_graph in the original)
#
/slurm/users/yangchentao/miniconda3/envs/verkko/bin/MBG \
  -i $iopt \
  -t $cpu \
  -k $kmer \
  -r 15000 -R 4000 \
  -w 100 \
  --kmer-abundance 1 \
  --unitig-abundance 2 \
  --error-masking=collapse-msat \
  --output-sequence-paths $outpre.paths.gaf \
  --out $outpre.hifi-resolved.gfa

#  Find coverage.
#    (rule hifi_coverage_csv in the original Snakefile)
#    (hifi-resolved.gfa -> hifi_nodecov.csv)
#    ($6 != "" is from 9e31a602925a477a7e52c277eda143e7bd20e52b)
#
awk 'BEGIN \
     { \
        FS="[ \t]+"; OFS="\t"; \
        print "node", "length", "coverage"; \
     } \
     $1 == "S" \
     { \
        if ($6 != "") {
          $4 = $6;
        }
        print $2, length($3), substr($4, 6); \
     }' \
< $outpre.hifi-resolved.gfa \
> $outpre.hifi_nodecov.csv
