#!/bin/sh
set -e

#  Figure out what reads we're going to use.  If correction is enabled, the
#  reads are in a single file that is a relative path away from us; but if it
#  is not enabled, the reads are in possibly multiple files that are at an
#  absolute path.
#

if [[ $# != 3   ]] ; then
	echo "Usage : $0 genome file_of_fastqs cpu"
	exit 1
fi

genome=$1
fof_fq=$2
cpu=$3

echo "Run juicer Hi-C pipeline...: " && date

wk=`pwd`
# build genome dir
if [ ! -d genome ]; then
	echo "genome dir not found, now build it"
	mkdir genome
	cd $wk/genome
	geno=`basename $genome`
	[ -e $geno ] || ln -s $genome ./genome.fasta
	cd $wk
else
	echo "genome dir found"
fi

# bwa index
if [ -e $wk/genome/genome.fasta.bwt ]; then
	echo "genome.fasta.bwt exists, skipping bwa index"
else
	bwa index $wk/genome/genome.fasta
fi

# restriction digest
if [ -e $wk/genome/genome.fasta__DpnII.txt ]; then
	echo "_DpnII.txt exist, skipping restriction digest"
else
	/usr/bin/python2.7 /slurm/users/yangchentao/software/T2T/juicer-v1.6/misc/generate_site_positions.py DpnII $wk/genome/genome.fasta $wk/genome/genome.fasta
fi

awk 'BEGIN{OFS = "\t"}{print $1,$NF}' $wk/genome/genome.fasta_DpnII.txt > $wk/genome/genome.fasta.sizes

# input fastq
if [ ! -d $wk/fastq ]; then
	mkdir fastq
	cat $fof_fq|while read a
	do
		ln -s $a ./fastq/
	done
else
	echo "fastq exists, do nothing"
fi


[ -d $wk/tmp ] || mkdir $wk/tmp

echo "#!/bin/bash
#SBATCH --job-name=juicer
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=200g
date
export JAVA_TOOL_OPTIONS=\"-Djava.io.tmpdir=$wk/tmp/ -XX:ParallelGCThreads=48\"
export TMPDIR=$wk/tmp/

bash /slurm/users/yangchentao/software/T2T/juicer-v1.6/scripts/juicer.sh -s DpnII -z $wk/genome/genome.fasta -y $wk/genome/genome.fasta_DpnII.txt -p $wk/genome/genome.fasta.sizes -D /slurm/users/yangchentao/software/T2T/juicer-v1.6/ -d $wk -t $cpu
date " > juicer.sh

echo "juicer preparation done, please run juicer.sh"

date
