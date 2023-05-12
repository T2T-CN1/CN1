#!/usr/bin/perl -w
use strict;

die "perl $0 <target.bed.complement> <filled_sequences.fasta.len>" unless @ARGV == 2;

my %hash;

open(IN, $ARGV[1]) or die $!;
while(<IN>){
	chomp;
	my ($id, $len) = (split /\t/)[0,1];
	$id =~ /(\w+):(\d+)-(\d+)/;
	my $start = $2;
	$hash{$start} = [$id, $len];
}
close IN;

my $n = keys %hash;
my $start = 1;
my $idx = 1;
open(IN, $ARGV[0]) or die $!;
while(<IN>){
	chomp;
	my ($scaf, $bg, $ed) = (split /\t/)[0,1,2];
	my $l = $ed-$bg;
	my $end = $start+$l-1;
	my $bg1 = $bg+1;
	print "chrY\t$start\t$end\t$idx\tW\t$scaf\t$bg1\t$ed\t+\n";
	$idx += 1;
	last if($idx == $n*2+2);
	$start = $end+1;
	my ($id, $len) = @{$hash{$ed}}[0,1];
	$end = $start+$len-1;
	print "chrY\t$start\t$end\t$idx\tW\t$id\t1\t$len\t+\n";
	$idx += 1;
	$start = $end+1;
}
close IN;

