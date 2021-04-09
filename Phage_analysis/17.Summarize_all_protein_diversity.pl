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

open OUT, ">All_phage_proteins/All_phage_protein_diversity_summary.txt";
print OUT "Plume ID\tAll phage protein num\tdRep95 all phage protein num\tAll phage protein cluster num\n";
foreach my $p_id (sort keys %PlumeID2){
	my $uknPhg = "All_phage_proteins/$p_id.all_phage_protein.faa";
	my $uknPhgdrep95 =  "All_phage_proteins/$p_id.all_phage.drep_rep_seq.fasta";
	my $uknPhgCluster =  "All_phage_proteins/$p_id.all_phage.drep95_looseCluster_rep_seq.fasta";
	my $uknPhg_SeqNum = `grep ">" $uknPhg | wc -l`; chomp $uknPhg_SeqNum; 
	my $uknPhgdrep95_SeqNum = `grep ">" $uknPhgdrep95 | wc -l`; chomp $uknPhgdrep95_SeqNum; 
	my $uknPhgCluster_SeqNum = `grep ">" $uknPhgCluster | wc -l`; chomp $uknPhgCluster_SeqNum;
	print OUT "$p_id\t$uknPhg_SeqNum\t$uknPhgdrep95_SeqNum\t$uknPhgCluster_SeqNum\n";
}
close OUT;