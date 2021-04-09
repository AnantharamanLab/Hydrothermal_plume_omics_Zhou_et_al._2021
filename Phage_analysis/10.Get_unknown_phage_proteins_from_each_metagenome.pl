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

my %MetaData = (); # Genome => [0] Taxonomy	[1] Environment	[2] Lysogenic property
open IN, "Gene_and_taxonomy_metadata.tsv";
while (<IN>){
	chomp;
	my @tmp = split (/\t/,$_);
	$MetaData{$tmp[0]}[0] = $tmp[1];
	$MetaData{$tmp[0]}[1] = $tmp[2];
	$MetaData{$tmp[0]}[2] = $tmp[3];
}
close IN;

`mkdir Unknown_phage_proteins`;

foreach my $p_id (sort keys %PlumeID){
	my %Seq = _store_seq("/slowdata/data1/hydrothermal_plume_omics/VIBRANT_running/VIBRANT_$p_id/VIBRANT_phages_$p_id/$p_id.phages_combined.faa");
	my %Seq2 = (); # Seq name changed
	foreach my $seq_name (sort keys %Seq){
		my $seq_name2 = $seq_name;
		$seq_name2 =~ s/ /\_/g;
		$seq_name2 =~ s/>/>$PlumeID{$p_id}~~/g;
		$Seq2{$seq_name2} = $Seq{$seq_name};
	}
	
	open OUT, ">Unknown_phage_proteins/$PlumeID{$p_id}.unknown_phage_protein.faa";
	foreach my $key (sort keys %Seq2){
		my ($key_gn) = $key =~ /^>(.+?)\_\d+?$/; #print "$MetaData{$key_gn}[0]\n";
		if ($MetaData{$key_gn}[0] eq "Unknown"){
			print OUT "$key\n$Seq2{$key}\n";
		}
	}
	close OUT;
}

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
