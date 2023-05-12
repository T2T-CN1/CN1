#!/usr/bin/perl -w

use strict;
use Getopt::Long;

unless(@ARGV){
	print "\n\tperl $0 <fasta file> < N size, default 1>\n\n";
	exit 0;
}


my $Nsize = $ARGV[1] || 1;

my $name = '';
my %fasta = ();

if($ARGV[0]=~/\.gz$/){
	open IN,"gunzip -c $ARGV[0] |" or die "$!\n";
}else{
	open IN,"$ARGV[0]" or die "$!\n";
}
while(my $line = <IN>){
	chomp $line;
	if($line =~ />(.*$)/){
		$name = (split /\s+/,$1)[0];
	}else{
		$fasta{$name} .= $line;
	}
}
close IN;


foreach my $name(sort keys %fasta){
	
	my $chrLen = length $fasta{$name};
	next if ($chrLen == 0);
	while($fasta{$name} =~ /N{$Nsize,}/gi){
		my $beg=$-[0]+1;
		my $end=$+[0];
		print "$name\t$beg\t$end\n";
	}
}
