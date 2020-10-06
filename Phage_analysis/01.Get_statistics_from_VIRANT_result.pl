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

my %Phage_sta = (); # plume id =>  s_id => len, lysogenic
my %MetaG_sta = (); # plume id =>  s_id => len
my %Phage_seq_num = (); # plume id =>  seq number 
my %MetaG_seq_num = (); # plume id =>  seq number
my %Phage_seq_nt_num = (); # plume id =>  seq number 
my %MetaG_seq_nt_num = (); # plume id =>  seq number
foreach my $plumeid (sort keys %PlumeID){
		my $phage_seq_num = 0;  my $metaG_seq_num = 0;
		my $phage_seq_nt_num = 0;  my $metaG_seq_nt_num = 0;
		my $phages_lysogenic_fna = "/slowdata/data1/hydrothermal_plume_omics/VIBRANT_running/VIBRANT_$plumeid/VIBRANT_phages_$plumeid/$plumeid.phages_lysogenic.fna";		
		my %Phages_lysogenic_seq = _store_seq($phages_lysogenic_fna);
		foreach my $s_id (sort keys %Phages_lysogenic_seq){
				my $len = length($Phages_lysogenic_seq{$s_id});
				$phage_seq_num++;  $phage_seq_nt_num += $len; 
				$Phage_sta{$PlumeID{$plumeid}}{$s_id} = "$len".","."lysogenic";
		}
		
		my $phages_lytic_fna = "/slowdata/data1/hydrothermal_plume_omics/VIBRANT_running/VIBRANT_$plumeid/VIBRANT_phages_$plumeid/$plumeid.phages_lytic.fna";		
		my %Phages_lytic_seq = _store_seq($phages_lytic_fna);
		foreach my $s_id (sort keys %Phages_lytic_seq){
				my $len = length($Phages_lytic_seq{$s_id});
				$phage_seq_num++;  $phage_seq_nt_num += $len; 
				$Phage_sta{$PlumeID{$plumeid}}{$s_id} = "$len".","."lytic";
		}
		
		my $metaG_fna = "/slowdata/data1/hydrothermal_plume_omics/VIBRANT_running/$plumeid.fasta";
		my %MetaG_seq = _store_seq($metaG_fna);
		foreach my $s_id (sort keys %MetaG_seq){
				my $len = length($MetaG_seq{$s_id});
				$metaG_seq_num++; $metaG_seq_nt_num += $len;
				$MetaG_sta{$PlumeID{$plumeid}}{$s_id} = $len;
		}		
		
		$Phage_seq_num{$PlumeID{$plumeid}} = $phage_seq_num;
		$MetaG_seq_num{$PlumeID{$plumeid}} = $metaG_seq_num;
		$Phage_seq_nt_num{$PlumeID{$plumeid}} = $phage_seq_nt_num;
		$MetaG_seq_nt_num{$PlumeID{$plumeid}} = $metaG_seq_nt_num;
		
}

open OUT, ">Phage_coding_density.tsv";
print OUT "Plume ID\tPhage sequence number\tMetagenome sequence number\tPhage sequence density (n / 1,000 metagenome sequences)\tPhage nucleotide number\tMetagenome nucleotide number\tPhage nucleotide density (n / 1,000,000 metagenome nucleotides)\n"; # Title
foreach my $p_id (sort keys %PlumeID2){
	my $density1 = 0; my $density2 = 0;
	#density1: sequence numbers
	$density1 = $Phage_seq_num{$p_id} / $MetaG_seq_num{$p_id};
	$density1 = 1000 * $density1; # normalized by 1,000 metagenome seq 
	$density1=sprintf "%.2f",$density1;
	#density2: sequence nt numbers
	$density2 = $Phage_seq_nt_num{$p_id} / $MetaG_seq_nt_num{$p_id};
	$density2 = 1000000 * $density2;# normalized by 1,000,000 metagenome nt
	$density2=sprintf "%.2f",$density2;
	print OUT "$p_id\t$Phage_seq_num{$p_id}\t$MetaG_seq_num{$p_id}\t$density1\t$Phage_seq_nt_num{$p_id}\t$MetaG_seq_nt_num{$p_id}\t$density2\n";
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
				($head) = $_ =~ /^(>.+?)\s/;
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
