#!/usr/bin/perl

use strict;
use warnings; 

#Separate mapping MetaG to bins that are from this environment.

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

my $MAG_adr = "/mnt/storage7/zhouzhichao/BinProject/hydrothermal_plume_omics/MAG_analysis";
my %CMD = (); # The command that cat all genome files (from one environment into one file); env => cmd
foreach my $env (@Env){
	$CMD{$env} = "cat ";
	foreach my $bin (sort keys %MAG_map) {
		if ($MAG_map{$bin}[1] eq $env){
			$CMD{$env} .= "$MAG_adr/\*_$bin\.genome ";
		}
	}
	$CMD{$env} .= "> $env.all.genome";
}

open OUT, ">tmp_All_genome_cat_metagenome_mapping.txt";
foreach my $env (@Env){
	`$CMD{$env}`;
	my $gn = "$env.all"; 
	system ("bowtie2-build $gn.genome $gn.genome_scaffold");
	
	foreach my $key (sort keys %map){
			if ($map{$key}[1] eq "DNA" && $map{$key}[4] eq $env){
				my $cmd = CMD($gn, $key);
				print OUT $cmd;
			}
	}
}
close OUT;
system ("cat tmp_All_genome_cat_metagenome_mapping.txt | /mnt/nfs_storage2/zhouzhichao/bin/parallel -j 24");
system ("rm tmp_All_genome_cat_metagenome_mapping.txt");
system ("rm *.all.genome_*bt2");

sub CMD{
my $v1 = $_[0];
my $v2 = $_[1];
my $cmd = "bowtie2 -x ${v1}.genome_scaffold  -1 /mnt/storage7/zhouzhichao/BinProject/QC_reads/${v2}_1.t.d.fastq  -2 /mnt/storage7/zhouzhichao/BinProject/QC_reads/${v2}_2.t.d.fastq  -S ${v1}.genome_$map{$v2}[0]_mapped.sam -p 128;";
$cmd = $cmd."samtools view -bS ${v1}.genome_$map{$v2}[0]_mapped.sam > ${v1}.genome_$map{$v2}[0]_mapped.bam -@ 128;";
$cmd = $cmd."samtools sort ${v1}.genome_$map{$v2}[0]_mapped.bam -o ${v1}.genome_$map{$v2}[0]_mapped.sorted.bam -@ 128;";
$cmd = $cmd."samtools index ${v1}.genome_$map{$v2}[0]_mapped.sorted.bam;";
$cmd = $cmd."samtools flagstat ${v1}.genome_$map{$v2}[0]_mapped.sorted.bam > ${v1}.genome_$map{$v2}[0]_mapped.sorted.stat;";
$cmd = $cmd."rm ${v1}.genome_$map{$v2}[0]_mapped.sam ${v1}.genome_$map{$v2}[0]_mapped.bam\n";
return $cmd;
}
