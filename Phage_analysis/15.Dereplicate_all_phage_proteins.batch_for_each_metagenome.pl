#!/usr/bin/perl

use strict;
use warnings;

my %PlumeID = (); # $p_id =>  Cayman.Deep
my %PlumeID2 = (); # Cayman.Deep => 1
open IN, "Plume_sample_id.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/,$_);
	$PlumeID{$tmp[0]} = $tmp[1];
	$PlumeID2{$tmp[1]} = 1;
}
close IN;

foreach my $p_id (sort keys %PlumeID2){
	my $input_phage_protein = "All_phage_proteins/$p_id.all_phage_protein.faa";
	`/slowdata/archive/anaconda2/bin/mmseqs easy-cluster $input_phage_protein All_phage_proteins/$p_id.all_phage.drep ./temp --min-seq-id 0.95 -c 0.8 -s 7.5 -e 0.001 --threads 10`;
	`rm -r temp`
	
	
}
