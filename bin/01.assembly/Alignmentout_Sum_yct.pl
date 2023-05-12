#!/usr/bin/perl 
#-------------------------------help-info-start--------------------------------#
=head1 Name 

  FILE: Alignmentout_Sum.pl

=head1 Description
 
   Obtain site information from multiple alignment file (in sam, bam, soap or map type).

=head1 Usage

  perl  Alignmentout_Sum.pl [options] -in <input file> -ref <input reference> -type <type>
  
  -in       <STR>   input alignment file
  -ref      <STR>   input reference sequence  (viral genome)
  -type     <STR>   type of input alignment file: sorted bam or pileup
  -outpre   [STR]   out prefix for results, default [myout]
  -outdir   [STR]   out directory for results, default [./]
  -dep      [INT]   only output site with depth more than int
  -heter    [boole] only output heterozygous sites
  -samtools [STR]   Samtools directory [defualt]
  -iTools   [STR]   iTools directory [defualt]
  -help     [STR]   print this help to screen

=head1 Example

  perl  Alignmentout_Sum.pl -in in.soap -ref ref.fa -type soap -filter yes -outpre myout -outdir ./

=head1 Contact

  Author: Chengran Zhou, zhouchengran@genomics.cn
          Jinmin Ma, majinmin@genomics.cn
  Organization: BGI-CNGB; BGI-SHENZHEN; SCU
  Version: 1.1

=cut
#-------------------------------help-info-end--------------------------------#

use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use lib "$Bin/";
use File::Basename qw(basename dirname);
use POSIX qw(strftime);

my ($Need_help,$in,$ref,$name_pre,$outdir,$type,$dep_cutoff,$only_heter,$samtools,$iTools);
GetOptions(
    "help!"      => \$Need_help,
	"in=s"		 => \$in,
	"ref=s"		 => \$ref,
	"outpre:s"   => \$name_pre,
	"outdir:s"   => \$outdir,
	"type=s"	 => \$type,
	"dep:i"      => \$dep_cutoff,
	"heter!"     => \$only_heter,
	"samtools:s"  => \$samtools,
	"iTools:s"   => \$iTools,
);

die `pod2text $0` if ($Need_help || !$in);
if ($type eq 'pileup'){
	die `pod2text $0` unless (defined $in);
}else{
	die `pod2text $0` unless (defined $in && defined $ref && defined $type);
}

#============================================================================#
#                              Global Variable                               #
#============================================================================#
$samtools ||= "samtools";
$iTools   ||= "iTools";

my $currentdir = `pwd`; chomp $currentdir;
$outdir ||= $currentdir;
mkdir $outdir unless (-d $outdir);

$name_pre ||= "myout";
$dep_cutoff ||= 3;
#my $pileup="$outdir/$name_pre.pileup";
my $pileup;

#============================================================================#
#                               Main process                                 #
#============================================================================#

system "$samtools faidx $ref" unless (-e "$ref.fai");
##format
if ($type =~/pileup/){
	$pileup=$in;
	
}elsif($type eq 'bam'){
	system "$samtools mpileup -AB -f $ref -d 10000 $in > $outdir/$name_pre.pileup  2>>$outdir/std.err";
	$pileup = "$outdir/$name_pre.pileup";	
}


##get results
open (IN, $pileup) or die $!;

##print results & entropy
open (OUT,">$outdir/$name_pre.base.xls") or die $!;
print OUT "Scaffold\tSite\tdepth\tdensity\tA\tT\tG\tC\tPi\tInsert\tDel\tReal depth\tEntropy\tmut_rate\tobserved depth\n";
#print OUT "Scaffold\tSite\tdepth\tdensity\tA\tT\tG\tC\tPi\n";
my %entropy;

# modify by yct 20200324
while (<IN>) {
	my (%base,$ndepth,$iddepth,$mut,$density, $npi);
	chomp;
	my ($NC,$site,$letter,$depth,$cov_base) = (split /\s+/, $_, 5);
	$letter = uc ($letter);
	if ($depth == 0){
		$cov_base="";
	}
	$cov_base =~ s/\^.|\$//g;
	while ($cov_base =~ /[+-]{1}(\d+)\w+/) {
		$cov_base =~ s/([+-]{1})$1(\w{$1})//;
		my $indel = uc($2);
		$base{$1}{$indel}++;
	}
	$base{A} = $cov_base =~ tr/Aa/Aa/;
	$base{T} = $cov_base =~ tr/Tt/Tt/;
	$base{G} = $cov_base =~ tr/Gg/Gg/;
	$base{C} = $cov_base =~ tr/Cc/Cc/;
	$base{$letter} = $cov_base =~ tr/,./,./;
	$base{skip}= $cov_base =~ tr/></></;
	my $nA=$base{A};
	my $nT=$base{T};
	my $nC=$base{C};
	my $nG=$base{G};
	$ndepth = $nA+$nT+$nC+$nG;
	##depth with indels
	$iddepth = $depth-$base{skip};
	my ($max_base) = (sort{$base{$b}<=>$base{$a}}  ('A','T','G','C',$letter));
	##mutation rate
	$mut = $iddepth ? sprintf "%0.4f",($iddepth - $base{$max_base})/($iddepth) : sprintf "%0.4f",0;
	##nucleotide diversity 	
	$npi = ($ndepth>1)? sprintf "%0.4f",(2*($nA*$nT+$nA*$nC+$nA*$nG+$nT*$nC+$nT*$nG+$nC*$nG)/(($ndepth)*($ndepth-1))) :sprintf "%0.4f",0;
	#$ndepth or $iddepth
	##density
	$density = ($ndepth>1)?sprintf "%0.4f",(log($ndepth)/log(2)):sprintf "%0.4f",0;

	my $Entropy = 0;
	my @tmp_output = ($NC, $site, $ndepth, $density);
	foreach my $base ('A','T','G','C') {
		push @tmp_output, $base{$base};
		next if ($iddepth==0);
		my $p = $base{$base}/$iddepth;
		next if $p == 0;
		$Entropy += $p * (log($p)/log(2));
	}
	push @tmp_output, $npi;
	# indel
	my%indel = ('+'=>'',
				'-'=>'',);
	foreach my $symbol ('+','-') {
		next unless exists $base{$symbol};
		next if ($iddepth==0);
		foreach my$indel (keys %{ $base{$symbol} }) {
			$indel{$symbol} .= "$indel($base{$symbol}{$indel});";
			my $p = $base{$symbol}{$indel}/$iddepth;
			$Entropy += $p * log($p)/log(2);
		}
		chop $indel{$symbol};
	}
	$Entropy *= -1;
	$Entropy =sprintf "%0.4f",$Entropy;
	push @tmp_output, ($indel{'+'}, $indel{'-'}, $iddepth, $Entropy, $mut, $depth);
	# only output heterozygous sites
	if ($only_heter){
		if ($ndepth >= $dep_cutoff && $Entropy > 0){
			print OUT join("\t", @tmp_output);
			print OUT "\n";
		}
	}else{
		if ($ndepth >= $dep_cutoff){
			print OUT join("\t", @tmp_output);
			print OUT "\n";
		}
	}
}

#system "rm $outdir/$name_pre.bam -f ";
if (-e "$outdir/$name_pre.sam") {
	system "rm $outdir/$name_pre.sam";
}
if (-e "$outdir/$name_pre.filter.sam"){
	system "rm $outdir/$name_pre.filter.sam";
}
if (-e "$outdir/$name_pre.filter.rmd.sam"){
	system "rm $outdir/$name_pre.filter.rmd.sam";
}
#system "rm $outdir/$name_pre.bam  $outdir/$name_pre.sam  $outdir/$name_pre.filter.sam  $outdir/$name_pre.filter.bam $outdir/$name_pre.filter.rmd.sam  -f";
close OUT;
