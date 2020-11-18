#!/usr/bin/perl

use strict;
use warnings; 

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
			$CMD{$env} .= "$MAG_adr/\*_$bin\.gene ";
		}
	}
	$CMD{$env} .= "> $env.all.gene";
}

open OUT, ">tmp_All_gene_cat_transcriptom_mapping.txt";
foreach my $env (@Env){
	`$CMD{$env}`;
	my $gn = "$env.all"; 
	system ("bowtie2-build $gn.gene $gn.gene_scaffold");
	

	foreach my $key (sort keys %map){
			if ($map{$key}[1] eq "cDNA" && $map{$key}[4] eq $env){
				if ($map{$key}[2] eq "SINGLE"){
					my $cmd = CMD1($gn, $key);
					print OUT $cmd;
				}elsif ($map{$key}[2] eq "PAIRED"){
					my $cmd = CMD2($gn, $key);
					print OUT $cmd;
				}
			}
	}
}
close OUT;
system ("cat tmp_All_gene_cat_transcriptom_mapping.txt | /mnt/nfs_storage2/zhouzhichao/bin/parallel -j 36");
system ("rm tmp_All_gene_cat_transcriptom_mapping.txt");
system ("rm *.gene_scaffold.*bt2");

#add parameters to let bowtie2 only find properly mapped reads (excluding single read mapping, unproperly mapping)
#and use the loose criterion to find mapping reads (which will report more potential reads)
#did not use the local mode, due to that I am afraid to find reads that is mapping to the boundary of genes and surpass the end of genes. (This point I am not sure about it)
#More information : Chinese explantation on parameters of bowtie2 : https://www.plob.org/article/4540.html
#English manual of details on the use of bowtie2 : http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#end-to-end-alignment-versus-local-alignment

#in the 1st test, I use "--no-mixed --no-discordant  --no-dovetail --no-contain --no-overlap --very-sensitive" to find mapped reads
#in the 2nd test, I use "--very-sensitive" to find mapped reads
sub CMD1{
my $v1 = $_[0];
my $v2 = $_[1];
my $cmd = "bowtie2 -x ${v1}.gene_scaffold  -U /mnt/storage7/zhouzhichao/BinProject/QC_reads/${v2}.trim_non_rRNA.fastq  -S ${v1}.gene_$map{$v2}[0]_mapped.sam -p 128  --very-sensitive;";
$cmd = $cmd."samtools view -bS ${v1}.gene_$map{$v2}[0]_mapped.sam > ${v1}.gene_$map{$v2}[0]_mapped.bam -@ 128;";
$cmd = $cmd."samtools sort ${v1}.gene_$map{$v2}[0]_mapped.bam -o ${v1}.gene_$map{$v2}[0]_mapped.sorted.bam -@ 128;";
$cmd = $cmd."samtools index ${v1}.gene_$map{$v2}[0]_mapped.sorted.bam;";
$cmd = $cmd."samtools flagstat ${v1}.gene_$map{$v2}[0]_mapped.sorted.bam > ${v1}.gene_$map{$v2}[0]_mapped.sorted.stat;";
$cmd = $cmd."rm ${v1}.gene_$map{$v2}[0]_mapped.sam ${v1}.gene_$map{$v2}[0]_mapped.bam\n";
return $cmd;
}

sub CMD2{
my $v1 = $_[0];
my $v2 = $_[1];
my $cmd = "bowtie2 -x ${v1}.gene_scaffold  -1 /mnt/storage7/zhouzhichao/BinProject/QC_reads/${v2}.trim_non_rRNA.R1.fastq  -2 /mnt/storage7/zhouzhichao/BinProject/QC_reads/${v2}.trim_non_rRNA.R2.fastq -S ${v1}.gene_$map{$v2}[0]_mapped.sam -p 128  --very-sensitive;";
$cmd = $cmd."samtools view -bS ${v1}.gene_$map{$v2}[0]_mapped.sam > ${v1}.gene_$map{$v2}[0]_mapped.bam -@ 128;";
$cmd = $cmd."samtools sort ${v1}.gene_$map{$v2}[0]_mapped.bam -o ${v1}.gene_$map{$v2}[0]_mapped.sorted.bam -@ 128;";
$cmd = $cmd."samtools index ${v1}.gene_$map{$v2}[0]_mapped.sorted.bam;";
$cmd = $cmd."samtools flagstat ${v1}.gene_$map{$v2}[0]_mapped.sorted.bam > ${v1}.gene_$map{$v2}[0]_mapped.sorted.stat;";
$cmd = $cmd."rm ${v1}.gene_$map{$v2}[0]_mapped.sam ${v1}.gene_$map{$v2}[0]_mapped.bam\n";
return $cmd;
}

