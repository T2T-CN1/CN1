#!/usr/bin/perl -w
use strict;

die "perl $0 <fasta>" unless @ARGV == 1;

open(IN, $ARGV[0]) or die $!;
$/ = ">"; <IN>;
while(<IN>){
	chomp;
	my ($id, $seq) = (split /\n/, $_, 2)[0,1];
	$seq =~ s/N//gi;
	$seq =~ s/\s+//g;
	print ">$id\n$seq\n";
}
close IN;

