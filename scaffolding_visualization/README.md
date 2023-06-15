### Scaffolding and visualization



## software preinstalled

- winnowmap2
- unimap
- LINKVIEW

```shell
# for hifiasm, maternal as an example
## assign to chromosome
xxx="your_data_path"
mkdir hifiasm
cd hifiasm
winnowmap -t 16 -W $xxx/00.dataset/t2t/chm13v2/repetitive_k19.txt -ax asm10 -H --MD $xxx/00.dataset/t2t/chm13v2/chm13v2.0.fasta  $xxx/CN1.hifiasm.dip.mat.p_ctg.fasta >CN1.hifiasm.mat.map2t2t.sam
paftools.js sam2paf -p CN1.hifiasm.mat.map2t2t.sam > CN1.hifiasm.mat.map2t2t.paf
python3 /bin/01.assembly/alignmentStatFromPaf.py CN1.hifiasm.mat.map2t2t.paf > CN1.hifiasm.mat.map2t2t.paf.stat
python3 /bin/01.assembly/assignUnitig2Chrs.py CN1.hifiasm.mat.map2t2t.paf.stat $xxx/CNhifiasm.dip.mat.p_ctg.fasta

## visualization alignment
cd splitbyChrs
# chr.list contains the list of chromsomes, one by line
cat  $xxx/chr.list | while read a
do
	[ -d $a ] || mkdir $a
	cd $a
	# generate karyotype file
	awk '$5=="+"' $a.hifiasm.unimap.paf|sort -k8n|cut -f 1,2,6,7|uniq > tmp
	python3 /bin/01.assembly/MakeLinkViewKaryotypeFromPaf.py tmp > karyotype.txt && rm tmp
	# linkview
	# make highlight.txt file, which contain the centromere and rDNA coordinates of CHM13
	python3 $software/Genome/LINKVIEW/LINKVIEW.py -t 3 -k karyotype.txt --svg_width 1800 --svg_height 600 --label_font_size 12 --label_angle 40 --chro_axis --gap_length 0.01 --svg2png_dpi 600 --no_dash  -hl highlight.txt $a.hifiasm.unimap.paf
	mv linkview_output.svg $a.hifiasm.svg
	mv linkview_output.png $a.hifiasm.png
	cd ..
done
```

