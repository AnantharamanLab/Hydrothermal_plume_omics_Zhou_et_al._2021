#!/usr/bin/perl

use strict;
use warnings; 

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

my $genome_dir = "/mnt/storage7/zhouzhichao/BinProject/hydrothermal_plume_omics/MAG_analysis";
my $cmd_cat = "cat "; # cat all Cayman bins
open IN, "Cayman_MAG.txt";
while (<IN>){
	chomp;
	$cmd_cat .= $genome_dir."/"."*_".$_.".genome ";
}
close IN;
$cmd_cat .= "> Cayman_genome_cat.genome";
system("$cmd_cat");


my $gn = "Cayman_genome_cat"; 
system ("bowtie2-build $gn.genome $gn.genome_scaffold");
open OUT, ">tmp_Cayman_genome_cat_metagenome_mapping.txt";
foreach my $key (sort keys %map){
        if ($map{$key}[1] eq "DNA" and $map{$key}[0] =~ /^Cym/){
			my $cmd = CMD($gn, $key);
			print OUT $cmd;
		}
}
close OUT;

system ("cat tmp_Cayman_genome_cat_metagenome_mapping.txt | parallel -j 50");
system ("rm tmp_Cayman_genome_cat_metagenome_mapping.txt");
system ("rm $gn.genome.*bt2");

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
