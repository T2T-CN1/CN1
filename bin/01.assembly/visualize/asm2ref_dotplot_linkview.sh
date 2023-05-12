#!/bin/sh

#!/usr/bin/bash
if [[ $# != 2   ]] ; then
	echo "Usage : $0 chr.list outpre"
	exit 1
fi

list=$1
outpre=$2

[ -d dotplot  ] || mkdir dotplot
[ -d linkview  ] || mkdir linkview

cat  chr.list | while read a
do
	[ -d $a  ] || mkdir $a
	cd $a
	sh /slurm/users/yangchentao//pipeline/customize/asm2ref/unimap/run_unimap_dotplot.sh /slurm/users/yangchentao/projects/T2T/00.dataset/t2t/chm13v2/chrs/$a.fa $a.fa $a.$outpre 48
	cp $a/$a.$outpre.dotplot.png dotplot
						
	# generate karyotype file
	awk '$5=="+"' $a.$outpre.unimap.paf|sort -k8n|cut -f 1,2,6,7|uniq > tmp
	python3 /slurm/users/yangchentao/projects/T2T/01.assembly/bin/MakeLinkViewKaryotypeFromPaf.py tmp > karyotype.txt && rm tmp
									
	# linkview
	## some chr has no rdna.bed.txt or chr1.unknown.highlight.txt, just ignore errors
	cat  cent.bed.txt rdna.bed.txt chr1.unknown.highlight.txt  > highlight.txt
	/slurm/users/yangchentao/miniconda3/bin/python3 /slurm/users/yangchentao/software/visualize/LINKVIEW/LINKVIEW.py -t 3 -k karyotype.txt --svg_width 1800 --svg_height 600 --label_font_size 12 --label_angle 40 --chro_axis --gap_length 0.01 --svg2png_dpi 600 --no_dash  -hl highlight.txt $a.$outpre.unimap.paf && echo "$a looks good!"
	mv linkview_output.svg $a.$outpre.svg
	mv linkview_output.png $a.$outpre.png
	cp $a/$a.$outpre.svg linkview
	cd ..
done

# archive
echo "archive plot results..."
tar czf dotplot.tgz dotplot
tar czf linkview.tgz linkview
echo "all done"
