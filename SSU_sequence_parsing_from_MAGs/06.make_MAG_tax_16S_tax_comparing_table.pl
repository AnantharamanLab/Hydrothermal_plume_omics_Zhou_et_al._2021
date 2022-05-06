#!/usr/bin/perl

use strict;
use warnings;

# store MAG tax
my %MAG_tax = (); # genome => GTDB curated tax  
open IN, "Hydrothermal_vent_ALR_GTDB-tk_result.txt";
while (<IN>){
	chomp;
	if (!/^Genome/){
		my @tmp = split (/\t/);
		$MAG_tax{$tmp[0]} = $tmp[7];
	}
}
close IN;

#Store bin ssu tax by Silva
my %SSu_tax = (); # scf_id => [0] lca_tax_embl_ebi_ena [1] lca_tax_gtdb [2] lca_tax_slv
open IN, "Bin_ssu_all.silva_result.txt";
while (<IN>){
	chomp;
	if (!/^sequence/){
		my @tmp = split (/\t/);
		$SSu_tax{$tmp[0]}[0] = $tmp[14];
		$SSu_tax{$tmp[0]}[1] = $tmp[15];
		$SSu_tax{$tmp[0]}[2] = $tmp[18];
	}
}
close IN;

#Store ssu scf info
my $head_ssu_info = "Scaffold ID\tHMM\ti-Evalue\tStart hit\tEnd hit\t16S/18S gene length\tRev. Complement\tSequence length\tGenome GTDB curated tax\tlca_tax_embl_ebi_ena\tlca_tax_gtdb\tlca_tax_slv";
my %SSu_info = (); # scaffold id => each line
open IN, "Total_scf_info.txt";
while (<IN>){
	chomp;
	if (!/^Sca/){
		my @tmp = split (/\t/);
		my $line_original = $_;
		my $scf_id = $tmp[0];
		my ($gn) = $scf_id =~ /^(.+?)\&\&/;
		
		my $line = $line_original;
		$line .= "\t".$MAG_tax{$gn}."\t".$SSu_tax{$scf_id}[0]."\t".$SSu_tax{$scf_id}[1]."\t".$SSu_tax{$scf_id}[2];
		$SSu_info{$scf_id} = $line;
	}
}
close IN;

# print out result
open OUT, ">MAG_tax_16S_tax_comparing_table.txt";
print OUT "$head_ssu_info\n";
foreach my $scf_id (sort keys %SSu_info){
	print OUT "$SSu_info{$scf_id}\n";
}
close OUT;