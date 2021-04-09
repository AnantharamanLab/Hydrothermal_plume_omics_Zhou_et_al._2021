#!/usr/bin/perl

use strict;
use warnings;

my %MetaData = (); # genome => [0] tax, [1] env, [2] lysogenic property [3] cluster

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

# Get taxonomy information
my %Tax_map = (); # Caudovirales_Ackermannviridae => Ackermannviridae
open IN, "Taxonomy_map.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/,$_);
	$Tax_map{$tmp[0]} = $tmp[1];
}
close IN;

# Get genome id hash
my %Genome_id = (); #
open IN, "/slowdata/data1/hydrothermal_plume_omics/VIBRANT_running/core-unique-flexible_result/core/genomes_list.mdf.txt";
while (<IN>){
	chomp;
	$Genome_id{$_} = 1;
}
close IN;

# Get hash for all sequences to cluster rep 
my %dRep95_2_Allseq = (); # dRep95 rep => protein seq \t join 
open IN, "All_plume_phage_drep95_cluster.tsv";
while (<IN>){
	chomp;
	my @tmp = split (/\t/,$_);
	$tmp[0] =~ s/Axial\.Plume\~\~/Axial\~\~/g; $tmp[0] =~ s/Axial\.Seawater\~\~/Axial\~\~/g;
	$tmp[1] =~ s/Axial\.Plume\~\~/Axial\~\~/g; $tmp[1] =~ s/Axial\.Seawater\~\~/Axial\~\~/g;
	if (!exists $dRep95_2_Allseq{$tmp[0]}){
		$dRep95_2_Allseq{$tmp[0]} = $tmp[1];
	}else{
		$dRep95_2_Allseq{$tmp[0]} .= "\t".$tmp[1];
	}
}
close IN;

=pod
my %Clr_rep_2_Allseq = (); # cluster rep = protein seq \t
open IN, "All_plume_phage_drep95_cluster.tsv";
while (<IN>){
	chomp;
	my @tmp = split (/\t/,$_); my $allseq = $dRep95_2_Allseq{$tmp[1]};
	if (!exists $Clr_rep_2_Allseq{$tmp[0]}){
		$Clr_rep_2_Allseq{$tmp[0]} = $allseq;
	}else{
		$Clr_rep_2_Allseq{$tmp[0]} .= "\t".$allseq;
	}
}	
close IN;
=cut

# Get protein to cluster hash
my %Gn2Clr21 = (); # Gn => Cluster => 1
my %Clr = (); # Cluster => 1
open IN, "ls /slowdata/data1/hydrothermal_plume_omics/VIBRANT_running/core-unique-flexible_result/core/*.faa |";
while (<IN>){
	chomp;
	my $file = $_;
	my ($clr) = $file =~ /\/slowdata\/data1\/hydrothermal_plume_omics\/VIBRANT_running\/core-unique-flexible_result\/core\/(.+?)\.faa/;
	$Clr{$clr} = 1;
	my %Seq = _store_seq1($file);
	foreach my $key (sort keys %Seq){
			my ($head) = $key =~ /^>(.+?)$/;
			my $tmp = "";
			if (!$dRep95_2_Allseq{$head}){
				print "$head\n";
			}else{
				$tmp = $dRep95_2_Allseq{$head}; 
			}
			my @All_seq = split (/\t/,$tmp); # get all the seq from the given cluster rep;
			foreach my $seq (@All_seq){
				my ($pro_gn) = $seq =~ /^(.+?)\_\d+?$/;	
				$Gn2Clr21{$pro_gn}{$clr} = 1;
			}			
	}
}
close IN;

# Get each genome cluster hits

foreach my $gn (sort keys %Genome_id){
	my $clr_hits_num = 0;
	foreach my $clr (sort keys %Clr){
		if ($Gn2Clr21{$gn}{$clr}){
			$clr_hits_num++;
		}
	}
	
	$MetaData{$gn}[3] = $clr_hits_num;
}

foreach my $p_id (sort keys %PlumeID){
	my $tax_file = "/slowdata/data1/hydrothermal_plume_omics/VIBRANT_running/VIVID_${p_id}_2020-10-05/$p_id.VIVID.virus-taxonomy.tsv";
	open IN, "$tax_file";
	while (<IN>){
		chomp;
		if (!/^scaffold/){
			my @tmp = split (/\t/,$_);
			my $genome_id = $tmp[0]; $genome_id =~ s/ /\_/g; $genome_id = $PlumeID{$p_id}.'~~'.$genome_id; $genome_id =~ s/Axial\.Plume\~\~/Axial\~\~/g; $genome_id =~ s/Axial\.Seawater\~\~/Axial\~\~/g;
			my $tmp = $tmp[2]."\_".$tmp[3]; my $tax = $Tax_map{$tmp};
			$MetaData{$genome_id}[0] = $tax;
			$MetaData{$genome_id}[1] = $PlumeID{$p_id}; # Store the metagenome origin a.k.a environment
		}
	}
	close IN;
	
	my $lysogenic_phage_file = "/slowdata/data1/hydrothermal_plume_omics/VIBRANT_running/VIBRANT_$p_id/VIBRANT_phages_$p_id/$p_id.phages_lysogenic.faa";
	my $lytic_phage_file = "/slowdata/data1/hydrothermal_plume_omics/VIBRANT_running/VIBRANT_$p_id/VIBRANT_phages_$p_id/$p_id.phages_lytic.faa";
	
	my %Lysogenic_phage_seq = _store_seq2($lysogenic_phage_file);
	my %Lytic_phage_seq = _store_seq2($lytic_phage_file);
	
	foreach my $head (sort keys %Lysogenic_phage_seq){
		my ($genome_id) = $head =~ /^>(.+?)\_\d+?$/; $genome_id =~ s/ /\_/g; $genome_id = $PlumeID{$p_id}.'~~'.$genome_id; $genome_id =~ s/Axial\.Plume\~\~/Axial\~\~/g; $genome_id =~ s/Axial\.Seawater\~\~/Axial\~\~/g;
		$MetaData{$genome_id}[2] = "lysogenic";
	}
	
	foreach my $head (sort keys %Lytic_phage_seq){
		my ($genome_id) = $head =~ /^>(.+?)\_\d+?$/; $genome_id =~ s/ /\_/g; $genome_id = $PlumeID{$p_id}.'~~'.$genome_id; $genome_id =~ s/Axial\.Plume\~\~/Axial\~\~/g; $genome_id =~ s/Axial\.Seawater\~\~/Axial\~\~/g;
		$MetaData{$genome_id}[2] = "lytic";
	}
}

# write Meta-data
open OUT, ">Gene_and_taxonomy_metadata.tsv";
print OUT "Genome\tTaxonomy\tEnvironment\tLysogenic property\tClusters\n";
foreach my $gn (sort keys %MetaData){
		print OUT "$gn\t$MetaData{$gn}[0]\t$MetaData{$gn}[1]\t$MetaData{$gn}[2]\t$MetaData{$gn}[3]\n";
}
close OUT;

sub _store_seq1{
	my $file = $_[0];
	my %Seq = (); my $head = "";
	open _IN, "$file";
	while (<_IN>){
		chomp;
		if (/>/){
			if (/\s/){
				($head) = $_ =~ /^(>.+?)$/;
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

sub _store_seq2{
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
