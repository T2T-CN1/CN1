#!/usr/bin/perl -w
use strict;
use FileHandle;

die "perl $0 <target.bed.seq> <outDir>" unless @ARGV == 2;

mkdir $ARGV[1] unless(-e $ARGV[1]);

my %fh;

open(IN, $ARGV[0]) or die $!;
while(<IN>){
	chomp;
	my ($scaf, $bg, $ed, $rid, $seq) = (split /\t/)[0,1,2,3,6];
	my $key = "$scaf-$bg-$ed";
	my $outDir = "$ARGV[1]/$key";
	mkdir $outDir unless(-e $outDir);
	
	if(!exists $fh{$key}){
		$fh{$key} = FileHandle->new("> $outDir/in.fasta");
	}
	if ($seq){
		$fh{$key} -> print(">$rid\n$seq\n");
	
	}else{
		print "[warning] $rid has no sequence, ignore...\n";
	}
}
close IN;

