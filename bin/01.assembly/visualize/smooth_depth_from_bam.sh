if [[ $# != 3  ]] ; then
	echo "Usage : $0 genome.fasta bam outprex"
	exit 1
fi

genome=$1
bam=$2
outpre=$3
win=1000

date
samtools faidx $genome
bedtools  makewindows -g $genome.fai -w 1000 >  $genome.fai.1k.bed
samtools depth -aa $bam |gzip >  $outpre.depth.gz 
zcat $outpre.depth.gz  | awk '{print $1"\t"$2-1"\t"$2"\t"$3}' | awk '$4>0'  | gzip >  $outpre.depth.bed.sort.gz
zcat $outpre.depth.bed.sort.gz | bedtools  map -a $genome.fai.1k.bed  -b - -c 4 -o median,mean,count >  $outpre.depth.1k.cov
date
