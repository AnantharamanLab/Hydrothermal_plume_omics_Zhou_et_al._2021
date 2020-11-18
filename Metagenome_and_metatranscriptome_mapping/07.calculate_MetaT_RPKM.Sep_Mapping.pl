#!/usr/bin/perl

use strict;
use warnings;

#1."per million" scaling factor is calculated as total number of mapped reads divided by 1,000,000
#2. Divide the read counts by the "per million" scaling factor. This normalizes for sequencing depth, giving you reads per million (RPM)
#3. Divide the RPM values by the length of the gene, in kilobases. This gives you RPKM.

my %map = (); # DNA and cDNA reads maps
open IN, "Transcriptom_map.mdf.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);  #$tmp[3] is the SRR id 
	$map{$tmp[3]}[0] = $tmp[0];   # SRR [0] => name
	$map{$tmp[3]}[1] = $tmp[1];   # SRR [1] => cDNA or DNA
	$map{$tmp[3]}[2] = $tmp[2];   # SRR [2] => PAIRED or SINGLE
	$map{$tmp[3]}[3] = $tmp[4];   # SRR [3] => reads numbers
	$map{$tmp[3]}[4] = $tmp[5];   # SRR [4] => Environment
}
close IN;

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


foreach my $env (@Env){
	my %h = (); # gene => MetaT name => RPKM
	my %MetaT = (); # Store the MetaT names

	open IN, "ls $env.all.gene*.C.*mapped.sorted.bam.pileup.out.parsed.csv | ";
	while (<IN>){
		chomp;
		my $file = $_ ;
		my ($name) = $_ =~ /\.all\.gene_(.+?)_mapped\.sorted\.bam\.pileup\.out\.parsed\.csv$/; # $env: the environment; $name: the name of MetaT library
		$MetaT{$name} =1;
		my $reads_num = 0;
		foreach my $key (sort keys %map){
			if ($map{$key}[0] eq $name){
				$reads_num = $map{$key}[3];
			}
		}
		
		open INN, $file;
		while (<INN>){
			chomp;
			my @tmp = split (/\,/,$_);
			my $gene_read_count = $tmp[1]; my $gene_read_len = $tmp[2]; 
			$h{$tmp[0]}{$name} = $gene_read_count * (1000000 / $reads_num) / ($gene_read_len /1000);
		}
		close INN;

	}
	close IN;

		
	#print table
	open OUT, ">$env.MetaT.RPKM.txt";
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
}