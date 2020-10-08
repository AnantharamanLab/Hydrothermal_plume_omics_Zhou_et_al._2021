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

my $all_phage_genomes = "";
foreach my $p_id (sort keys %PlumeID){
	my %Seq = _store_seq("/slowdata/data1/hydrothermal_plume_omics/VIBRANT_running/VIBRANT_$p_id/VIBRANT_phages_$p_id/$p_id.phages_combined.faa");
	my %Seq2 = (); # Seq name changed
	foreach my $seq_name (sort keys %Seq){
		my $seq_name2 = $seq_name;
		$seq_name2 =~ s/ /\_/g;
		$seq_name2 =~ s/>/>$PlumeID{$p_id}~~/g;
		$Seq2{$seq_name2} = $Seq{$seq_name};
	}
	open OUT, ">/slowdata/data1/hydrothermal_plume_omics/VIBRANT_running/VIBRANT_$p_id/VIBRANT_phages_$p_id/$p_id.phages_combined.mdf.faa";
	foreach my $key (sort keys %Seq2){
			print OUT "$key\n$Seq2{$key}\n";
	}
	close OUT;
	$all_phage_genomes .= "/slowdata/data1/hydrothermal_plume_omics/VIBRANT_running/VIBRANT_$p_id/VIBRANT_phages_$p_id/$p_id.phages_combined.mdf.faa ";
}

`cat $all_phage_genomes > All_plume_phage_proteins.faa`;

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