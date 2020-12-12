#!/usr/bin/perl

use strict;
use warnings;

my %MAG_map = (); # bin => [0] Microbial group; [1] Genome Phylogeny; [2] Strain name

open IN,"MAG_info.txt";
while (<IN>){
	chomp;
	if (!/^#/){
		my @tmp = split (/\t/,$_);
		$MAG_map{$tmp[0]}[0] = $tmp[2];
		$MAG_map{$tmp[0]}[1] = $tmp[4];
		$MAG_map{$tmp[0]}[2] = $tmp[5];
	}
}
close IN;

open IN, "Bin.list";
while (<IN>){
	chomp;
	my $bin = $_;
	print "$bin\t$MAG_map{$bin}[0]\t$MAG_map{$bin}[1]\t$MAG_map{$bin}[2]\n";
}
close IN;