#!/bin/sh
[ -d ragtag_linkview ] || mkdir ragtag_linkview
[ -d ragtag_dotplot  ] || mkdir ragtag_dotplot
#!/usr/bin/bash
if [[ $# != 2   ]] ; then
	echo "Usage : $0 chr.list outpre"
	exit 1
fi

list=$1
outpre=$2

 cat  $list | while read a
do
	[ -d $a  ] || mkdir $a
	cd $a/ragtag
	if [ -e ragtag.patch.fasta ]; then
		echo "$a ragtag succeed"

		sh /slurm/users/yangchentao/pipeline/customize/asm2ref/unimap/run_unimap_dotplot.sh /slurm/users/yangchentao/projects/T2T/00.dataset/t2t/chm13v2/chrs/$a.fa ragtag.patch.fasta $a.$outpre 48
		cp $a.$outpre.dotplot.png ../../ragtag_dotplot
						
		# generate karyotype file
		cat $a.$outpre.unimap.paf|sort -k8n|cut -f 1,2,6,7|uniq > tmp
		python3 /slurm/users/yangchentao/projects/T2T/01.assembly/bin/MakeLinkViewKaryotypeFromPaf.py tmp > karyotype.txt && rm tmp
									
		# linkview
		## some chr has no rdna.bed.txt or chr1.unknown.highlight.txt, just ignore errors
		cp /slurm/users/yangchentao/projects/T2T/01.assembly/F1/FF11/assembly_versions/draft_v0.0/maternal/verkko/splitbyChrs/$a/highlight.txt .
		/slurm/users/yangchentao/miniconda3/bin/python3 /slurm/users/yangchentao/software/visualize/LINKVIEW/LINKVIEW.py -t 3 -k karyotype.txt --svg_width 1800 --svg_height 600 --label_font_size 12 --label_angle 40 --chro_axis --gap_length 0.01 --svg2png_dpi 600 --no_dash  -hl highlight.txt $a.$outpre.unimap.paf && echo "$a looks good!"
		mv linkview_output.svg $a.$outpre.svg
		mv linkview_output.png $a.$outpre.png
		cp $a.$outpre.svg ../../ragtag_linkview
	else
		echo "$a ragtag failed"
	fi
	cd ../..
done

# archive
echo "archive plot results..."
tar czf ragtag_dotplot.tgz ragtag_dotplot
tar czf ragtag_linkview.tgz ragtag_linkview
echo "all done"
