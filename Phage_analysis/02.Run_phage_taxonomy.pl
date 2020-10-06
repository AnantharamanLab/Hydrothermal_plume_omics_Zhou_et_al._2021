#!/usr/bin/perl

use strict;
use warnings;

my %PlumeID = (); # Cayman_Deep_k45_95_10_min1000 =>  Cayman.Deep
my %PlumeID2 = (); # Cayman.Deep => 1
open IN, "Plume_sample_id.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/,$_);
	$PlumeID{$tmp[0]} = $tmp[1];
	$PlumeID2{$tmp[1]} = 1;
}
close IN;

foreach my $p_id (sort keys %PlumeID){
	`python3 /slowdata/data4/VIVID/phage_taxonomy_tool_v6.py -i /slowdata/data1/hydrothermal_plume_omics/VIBRANT_running/VIBRANT_${p_id}/VIBRANT_phages_${p_id}/${p_id}.phages_combined.faa -t 10`;
}


