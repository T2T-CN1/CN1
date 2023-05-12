#!/usr/bin/perl -w
use strict;

die "perl $0 <target.bed> <inDir>" unless @ARGV == 2;

open(IN, $ARGV[0]) or die $!;
while(<IN>){
	chomp;
	my ($scaf, $bg, $ed) = (split /\t/)[0,1,2];
	my $key = "$scaf-$bg-$ed";
	my $path = "$ARGV[1]/$key/in.fasta.out.cons.removeN";
	if(-e $path and -s $path){
		open(IN1, $path);
		$/ = ">"; <IN1>;
		while(<IN1>){
			my ($id, $seq) = (split /\n/, $_, 2)[0,1];
			$seq =~ s/\s+//g;
			print ">$scaf:$bg-$ed\n$seq\n";
		}
		close IN1;
		$/ = "\n";
	}
	else{
		print STDERR "$key/in.fasta.out.cons.removeN not found\n";
	}
}
close IN;

