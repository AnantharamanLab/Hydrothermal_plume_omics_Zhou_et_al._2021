#!/usr/bin/perl

use strict;
use warnings;


#TPM (Transcripts Per Kilobase Million)
#1. Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
#2. Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
#3. Divide the RPK values by the “per million” scaling factor. This gives you TPM.

my %map = (); # DNA and cDNA reads maps
open IN, "Transcriptom_map.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	$map{$tmp[3]}[0] = $tmp[0];   # SRR [0] => name
	$map{$tmp[3]}[1] = $tmp[1];   # SRR [1] => cDNA or DNA
	$map{$tmp[3]}[2] = $tmp[2];   # SRR [2] => PAIRED or SINGLE
	$map{$tmp[3]}[3] = $tmp[4];   # SRR [3] => reads numbers
}
close IN;

my %h = (); my %MetaT = ();
open IN, "ls Cayman_gene_cat.gene_*_mapped.sorted.bam.pileup.out.parsed.csv | ";
while (<IN>){
	chomp;
	my $file = $_;
	my ($name) = $_ =~ /^Cayman\_gene\_cat\.gene\_(.+?)\_mapped\.sorted\.bam\.pileup\.out\.parsed\.csv$/;
	$MetaT{$name} =1;

	my $total_rpk = 0;
	my %Gene_ids = (); # Store the gene id in Cayman_gene_cat.gene_*_mapped.sorted.bam.pileup.out.parsed.csv
	open INN, $file;
	while (<INN>){
		chomp;
		my @tmp = split (/\,/,$_);
		my $gene_read_count = $tmp[1]; my $gene_read_len = $tmp[2]; 
		my $rpk = $gene_read_count / ($gene_read_len /1000);
		$total_rpk += $rpk;
		$h{$tmp[0]}{$name} = $rpk;
		$Gene_ids{$tmp[0]} = 1;
	}
	close INN;
	
	my $scaling_factor = $total_rpk / 1000000;
	
	foreach my $gene (sort keys %Gene_ids){
		my $tmp = $h{$gene}{$name};
		if ($tmp){
			$h{$gene}{$name} = $tmp / $scaling_factor;
		}
	}

}
close IN;
	
#print table
open OUT, ">MetaT.TPM.txt";
my $row=join("\t", sort keys %MetaT);
print OUT "Head\t$row\n";
foreach my $tmp1 (sort keys %h)
{
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (sort keys %MetaT)
        {
                if (exists $h{$tmp1}{$tmp2})
                {
                        push @tmp, $h{$tmp1}{$tmp2};
                }
                else
                {
                        push @tmp,"0"
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;
