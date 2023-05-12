if [[ $# < 2  ]] ; then
	echo "Usage : $0 outpre assembly.fa [ref.fa]"
	exit 1
elif [[ $# == 3 ]]; then
	outpre=$1
	que=$2
	ref=$3
	echo "date
python3 /hwfsxx1/ST_HN/P18Z10200N0197/yangchentao/software/T2T/quast-5.1/quast.py -o  ${outpre}_quast -r $ref -t 24 $que
date " > ${outpre}_quast.sh

	sel-qsub evo 40g 24 ${outpre}_quast.sh
elif [[ $# == 2 ]]; then
	outpre=$1
	que=$2
	echo "date
python3 /hwfsxx1/ST_HN/P18Z10200N0197/yangchentao/software/T2T/quast-5.1/quast.py -o  ${outpre}_quast  -t 24 $que
date " > ${outpre}_quast.sh

	sel-qsub evo 40g 24 ${outpre}_quast.sh
fi
