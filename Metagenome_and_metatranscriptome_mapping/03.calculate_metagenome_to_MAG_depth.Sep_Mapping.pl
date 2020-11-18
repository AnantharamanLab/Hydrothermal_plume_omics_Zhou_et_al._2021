#!/usr/bin/perl

use strict;
use warnings;

my %MAG_map = (); # bin => [0] Microbial group; [1] Environment
my @Env = ();my %Env = (); # Store the Environments; 
open IN,"MAG_info.txt";
while (<IN>){
	chomp;
	if (!/^#/){
		my @tmp = split (/\t/,$_);
		$MAG_map{$tmp[0]}[0] = $tmp[2];
		$MAG_map{$tmp[0]}[1] = $tmp[-1];
		if (!exists $Env{$tmp[-1]}){
			$Env{$tmp[-1]} = 1; push @Env,$tmp[-1];
		}
	}
}
close IN;

#  $env.all.genome_LBKM.D.NBP_mapped.sorted.bam

my %CMD = (); # The command that calculate all depth files by jgi_summarize_bam_contig_depths; env => cmd
foreach my $env (@Env){
	$CMD{$env} = "jgi_summarize_bam_contig_depths --outputDepth $env.all.genome.depth.txt  --pairedContigs $env.all.genome.paired.txt ";
}
open IN, "ls *.D.*mapped.sorted.bam|";
while (<IN>){
	chomp;
	my $bam_file  = $_;
	my ($env1) = $_ =~ /^(.+?)\.all\.genome_/;
	foreach my $env2 (sort keys %CMD){
		if ($env1 eq $env2){
			$CMD{$env2} .= $bam_file." ";
		}
	}
}
close IN;

#test output
open OUT, ">tmp.calculate_metagenome_to_MAG_depth.Sep_Mapping.sh";
foreach my $env (sort keys %CMD){
	print OUT "$CMD{$env}\n";
}
close OUT;

system ("cat tmp.calculate_metagenome_to_MAG_depth.Sep_Mapping.sh | /mnt/nfs_storage2/zhouzhichao/bin/parallel -j 10");
system ("rm tmp.calculate_metagenome_to_MAG_depth.Sep_Mapping.sh");


