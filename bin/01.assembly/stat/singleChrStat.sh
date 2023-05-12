#!/bin/sh

echo -e "chr\t>1Mb\t<1Mb\tcompleteness"
cat  /slurm/users/yangchentao/projects/T2T/01.assembly/chr.list | while read a 
do
	[ -d $a ] || mkdir $a
	cd $a
	if [ -e karyotype.txt ] ; then
		large=`head -1 karyotype.txt|xargs -n1|awk -F ":" '$3>=1000000' |wc -l`
		small=`head -1 karyotype.txt|xargs -n1|awk -F ":" '$3<1000000' |wc -l`
		comp=`python3 /slurm/users/yangchentao/projects/T2T/01.assembly/bin/RefCoverageStatFromPaf.py $a.salsa.unimap.paf|cut -f 4`
		echo -e "$a\t$large\t$small\t$comp"
	fi

	cd ..
done

