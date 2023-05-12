#!/usr/bin/bash
if [[ $# != 2  ]] ; then
	echo "Usage : $0 input_bam outdir"
	exit 1
fi

fbam=$1
outdir=$2

if [ ! -e ${fbam}.bai ]; then
	echo "can not find index file of $fbam"
	exit 1
fi

[ -d $outdir ] || mkdir $outdir
[ -d tmp_shell ] || mkdir tmp_shell

wk=`pwd`
samtools view -H $fbam |grep '^@SQ' | cut -f 2,3 |sed 's/SN://g;s/LN://g' > $wk/tmp_shell/job_resources.txt

# make shell
for i in `samtools view -H $fbam|grep '^@SQ' |cut -f2|sed 's/SN://g'`
do
	echo "date
samtools view -F 256 -h -bS -o $wk/$outdir/$i.bam $wk/$fbam $i  
samtools index $wk/$outdir/$i.bam
python3 /hwfsxx1/ST_HN/P18Z10200N0197/yangchentao/software/T2T/NucFreq/NucPlot.py -y 200 -t 8 $wk/$outdir/$i.bam $wk/$outdir/$i.nucplot.png
date " > $wk/tmp_shell/$i.sh
done

# submit jobs
cd $wk/tmp_shell
echo "start qsub jobs...."
cat job_resources.txt|while read scaf length
do
	if [ $length -gt 100000000 ]; then
		sel-qsub evo 50g 8 $scaf.sh
	elif [ $length -gt 50000000 ]; then
		sel-qsub evo 15g 8 $scaf.sh
	elif [ $length -gt 10000000 ]; then
		sel-qsub evo 5g 8 $scaf.sh
	else
		sel-qsub evo 3g 8 $scaf.sh
	fi
done

cd $wk
