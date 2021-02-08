#!/usr/bin/perl

use strict;
use warnings;

my %PlumeID = (); # Cayman_Deep_k45_95_10_min1000 => Cayman.Deep
my %PlumeID2 = (); # Cayman.Deep => 1
open IN, "Plume_sample_id.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/,$_);
	$PlumeID{$tmp[0]} = $tmp[1];
	$PlumeID2{$tmp[1]} = 1;
}
close IN;

my %Genome_list = (); # genome id => 1
foreach my $p_id (sort keys %PlumeID){
	my %Seq = _store_seq("/slowdata/data1/hydrothermal_plume_omics/VIBRANT_running/VIBRANT_$p_id/VIBRANT_phages_$p_id/$p_id.phages_combined.mdf.faa");
	foreach my $head (sort keys %Seq){
			my ($gn) = $head =~ /^>(.+?)\_\d+?$/;
			$Genome_list{$gn} = 1;
	}
}

open OUT, ">genomes_list.txt";
foreach my $key (sort keys %Genome_list){
		print OUT "$key\n";
}
close OUT;

sub _store_seq{
	my $file = $_[0];
	my %Seq = (); my $head = "";
	open _IN, "$file";
	while (<_IN>){
		chomp;
		if (/>/){
			if (/\s/){
				($head) = $_ =~ /^(>.+?)\t/;
				$Seq{$head} = "";
			}else{
				($head) = $_ =~ /^(>.+?)$/;
				$Seq{$head} = "";
			}
		}else{
			$Seq{$head} .= $_;
		}
	}
	close _IN;
	return %Seq;
}