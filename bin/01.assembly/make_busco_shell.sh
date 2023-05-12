if [[ $# != 2 ]] ; then
	echo "Usage : $0 assembly.fa outpre"
	exit 1
fi

asm=$1
outpre=$2
db='/hwfsxx1/ST_HN/P18Z10200N0197/yangchentao/T2T/00.dataset/busco_primates_odb10/'
cpu=24

echo "date
source ~/.bash_condainit
conda activate busco
busco -i $asm -o ${outpre}_busco -l $db -m geno --cpu $cpu --offline && echo 'done'
date " > ${outpre}_busco.sh
sel-qsub evo 50g 24 ${outpre}_busco.sh
