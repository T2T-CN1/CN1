# Telomere-to-telomere of CN1 Han Chinese



**INPORTANT NOTE**

All datasets are prepared and ready to release, however, becasuse of policy issue, we need to hold 1-2 month for permission then to make it piblic. 

All raw sequence data used in this study are available for CN1 at GitHub: https://github.com/T2T-CN1/CN1. The reads are also deposited in the CNCB under accession number HRA004405. The final CN1 curated assemblies are available in the CNCB with the accession number GWHCBHP00000000 under the BioProject ID PRJCA016397. The accession numbers for the maternal and paternal haplotypes are GWHCBHM00000000 and GWHCBHQ00000000 in CNCB. We also submitted the assemblies to CNGB with the accession numbers of CNA0069006-CNA0069008 for combined, maternal, and paternal genome, respectively. The sample datasets are available on CNGB with the BioProject ID CNP0004252. Assemblies, annotation, and variant results are available at https://genome.zju.edu.cn.



## Customized codes 

- scaffolding and visualization

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

  

- Heterozygosity analysis

  ```shell
  # MF2v0.8. maske cent
  mat_fa=xxx/v0.8/MF2_mat.v0.8.fasta
  pat_fa=xxx/v0.8/MF2_pat.v0.8.fasta
  
  # split to chr
  fastaKit -ap mat -o MF2_mat.v0.8.addTag.fasta $mat_fa
  fastaKit -ap pat -o MF2_pat.v0.8.addTag.fasta $pat_fa
  python3 bin/split_genome_bychrs.py MF2_mat.v0.8.addTag.fasta MF2_pat.v0.8.addTag.fasta
  
  # autosome
  for i in `seq 1 22`
  do
   chr="chr"$i
   [ -d $chr ] || mkdir $chr
   cd $chr
   echo "#!/bin/bash
  date
  minimap2 -x asm5  -t 24 --cs mat_${chr}.fa pat_${chr}.fa |sort -k6,6 -k8,8n  > $chr.minimap.sort.paf
  paftools.js call $chr.minimap.sort.paf > $chr.var.txt
  date" >  $chr.minimap.sh
   #sbatch -c 24 --mem 16g $chr.minimap.sh
   grep -v '^R' $chr.var.txt |cut -f 2-4|sed 's/mat_//g' > $chr.var.bed
   sort -k1,1V -k2n $chr.var.bed |bedtools merge -i - > $chr.var.bed.merge
   cd ..
  done
  
  cut -f 1,2 xxx/v0.8/MF2_mat.v0.8.fasta.fai > mat.v0.8.genome.size
  bedtools makewindows -g mat.v0.8.genome.size  -w 500000 -s 500000 > mat.genome.500k.bed
  cat chr*/*.var.txt |grep -v '^R' >autosome.allvar.txt
  cat chr*/*.var.bed.merge > autosome.allvar.bed
  bedtools coverage -a ../2.collect/mat.genome.500k.bed -b autosome.allvar.bed > autosome.500k.allvar.cov
  python3 /bin/makeGFA_from_Bed.py autosome.500k.allvar.cov 0.0004 > autosome.500k.allvar.0.0004.ab.gfa
  
  # then visualize autosome.500k.allvar.0.0004.ab.gfa using Bandage
  ```

  

- SV manual check

  [adopted version of bamsnap](https://github.com/zy041225/bamsnap)

  
