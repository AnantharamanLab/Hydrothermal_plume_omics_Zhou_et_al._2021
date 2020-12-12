#!/usr/bin/perl

use strict;
use warnings;

=pod
-----------------------------------------
hmm	131
custom_blastp	8
KEGG_function	26
dbCAN	145
MEROPS	172
-----------------------------------------
=cut

#hmm dir; result files are ended with .hmmscan.result.txt
my $hmm_dir = "/mnt/storage7/zhouzhichao/BinProject/hydrothermal_plume_omics/MAG_analysis/do_hmm_annotation";
#custom_blastp dir; result files are ended with .custom_blastp.result.txt
my $blastp_dir = "/mnt/storage7/zhouzhichao/BinProject/hydrothermal_plume_omics/MAG_analysis/custom_blastp";
#KEGG_function dir; result files are ended with .KEGG_function.result.txt
my $KEGG_dir = "/mnt/storage7/zhouzhichao/BinProject/hydrothermal_plume_omics/MAG_analysis/parse_KEGG_result";
#dbCAN dir; result files are ended with .dbCAN.result.txt
my $dbCAN_dir = "/mnt/storage7/zhouzhichao/BinProject/hydrothermal_plume_omics/MAG_analysis/dbCAN";
#MEROPS dir; result files are ended with .MEROPS.result.txt
my $MEROPS_dir = "/mnt/storage7/zhouzhichao/BinProject/hydrothermal_plume_omics/MAG_analysis/MEROPS/Extracellular_analysis";

my %MetaT_RPKM = (); # gene => MetaT => RPKM

my @head = ();
open IN, "/mnt/storage7/zhouzhichao/BinProject/hydrothermal_plume_omics/Cayman_Mapping.MetaT_new/MetaT.RPKM.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/,$_);
	if ($tmp[0] eq "Head"){
		@head = @tmp;
	}else{
		my $gene = $tmp[0];  
		for(my $i=1; $i<=$#head; $i++){
			$MetaT_RPKM{$gene}{$head[$i]} = $tmp[$i];
		}
	}	
}
close IN;

my %MetaT_id = (); my @MetaT_id = (); 
open IN,"MetaT_list.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/,$_);
	if ($tmp[0] =~ /^Cym/){
		$MetaT_id{$tmp[0]} = 1;
		push @MetaT_id, $tmp[0];	
	}
}
close IN;

#-----------------------------------
my %gene_id = (); # Store all genes
my %Fun_id = (); # Store the function ids
my @Fun_id = ();
my %Bin_id = (); # Store the bin (aka. genome) names
my %Fun_result = (); # fun => bin => genes
open IN, "ls $hmm_dir/*SZU*.hmmscan.result.txt | ";
while (<IN>){
	chomp;
	my $file = $_;
	my ($bin) = $_ =~ /.+?\_(SZU.+?)\.hmmscan\.result\.txt/; $Bin_id{$bin} = 1;
	open INN, $file;
	while (<INN>){
		chomp;
		my @tmp = split (/\t/);
		my $fun = $tmp[0]; 
		if (!exists $Fun_id{$fun}){
			$Fun_id{$fun} = 1;
			push @Fun_id, $fun;
		}
		my @genes = ();
		for(my $i=1; $i<=$#tmp; $i++){
			push @genes, $tmp[$i];
			if ($tmp[$i]){
				$gene_id{$tmp[$i]} = 1;
			}
		}
		$Fun_result{$fun}{$bin} = join("\t",@genes);
	}
	close INN;
}
close IN;

open IN, "ls $blastp_dir/*SZU*.custom_blastp.result.txt | ";
while (<IN>){
	chomp;
	my $file = $_;
	my ($bin) = $_ =~ /.+?\_(SZU.+?)\.custom_blastp\.result\.txt/; 
	open INN, $file;
	while (<INN>){
		chomp;
		my @tmp = split (/\t/);
		my $fun = $tmp[0]; 
		if (!exists $Fun_id{$fun}){
			$Fun_id{$fun} = 1;
			push @Fun_id, $fun;
		}
		my @genes = ();
		for(my $i=1; $i<=$#tmp; $i++){
			push @genes, $tmp[$i];
			if ($tmp[$i]){
				$gene_id{$tmp[$i]} = 1;
			}			
		}
		$Fun_result{$fun}{$bin} = join("\t",@genes);
	}
	close INN;
}
close IN;

open IN, "ls $KEGG_dir/*SZU*.KEGG_function.result.txt | ";
while (<IN>){
	chomp;
	my $file = $_;
	my ($bin) = $_ =~ /.+?\_(SZU.+?)\.KEGG_function\.result\.txt/; 
	open INN, $file;
	while (<INN>){
		chomp;
		my @tmp = split (/\t/);
		my $fun = $tmp[0]; 
		if (!exists $Fun_id{$fun}){
			$Fun_id{$fun} = 1;
			push @Fun_id, $fun;
		}
		my @genes = ();
		for(my $i=1; $i<=$#tmp; $i++){
			push @genes, $tmp[$i];
			if ($tmp[$i]){
				$gene_id{$tmp[$i]} = 1;
			}			
		}
		$Fun_result{$fun}{$bin} = join("\t",@genes);
	}
	close INN;
}
close IN;

open IN, "ls $dbCAN_dir/*SZU*.dbCAN.result.txt | ";
while (<IN>){
	chomp;
	my $file = $_;
	my ($bin) = $_ =~ /.+?\_(SZU.+?)\.dbCAN\.result\.txt/; 
	open INN, $file;
	while (<INN>){
		chomp;
		my @tmp = split (/\t/);
		my $fun = $tmp[0]; 
		if (!exists $Fun_id{$fun}){
			$Fun_id{$fun} = 1;
			push @Fun_id, $fun;
		}
		my @genes = ();
		for(my $i=1; $i<=$#tmp; $i++){
			push @genes, $tmp[$i];
			if ($tmp[$i]){
				$gene_id{$tmp[$i]} = 1;
			}			
		}
		$Fun_result{$fun}{$bin} = join("\t",@genes);
	}
	close INN;
}
close IN;

open IN, "ls $MEROPS_dir/*SZU*.MEROPS.result.txt | ";
while (<IN>){
	chomp;
	my $file = $_;
	my ($bin) = $_ =~ /.+?\_(SZU.+?)\.MEROPS\.result\.txt/; 
	open INN, $file;
	while (<INN>){
		chomp;
		my @tmp = split (/\t/);
		my $fun = $tmp[0]; 
		if (!exists $Fun_id{$fun}){
			$Fun_id{$fun} = 1;
			push @Fun_id, $fun;
		}
		my @genes = ();
		for(my $i=1; $i<=$#tmp; $i++){
			push @genes, $tmp[$i];
			if ($tmp[$i]){
				$gene_id{$tmp[$i]} = 1;
			}			
		}
		$Fun_result{$fun}{$bin} = join("\t",@genes);
	}
	close INN;
}
close IN;

#-----------------------------------

#-----------------------------------





my %Fun2MetaT = (); #Fun => MetaT => added in total MetaT RPKM
foreach my $fun (sort keys %Fun_id){
	foreach my $MetaT (sort keys %MetaT_id){
		my $MetaT_RPKM_value = 0;		
		foreach my $bin (sort keys %Bin_id){
			foreach my $gene (sort keys %gene_id){
				if ($Fun_result{$fun}{$bin} && $Fun_result{$fun}{$bin} =~ /$gene/){   # my %Fun_result = (); # fun => bin => genes
					if ($MetaT_RPKM{$gene}{$MetaT}){
						$MetaT_RPKM_value += $MetaT_RPKM{$gene}{$MetaT}; # my %MetaT_RPKM = (); # gene => MetaT => RPKM
					}else{
						$MetaT_RPKM_value += 0;
					}					
				}
			}		
		}	
		$Fun2MetaT{$fun}{$MetaT} = $MetaT_RPKM_value;
	}
}


open OUT, ">Fun2MetaT_abundance.txt";
my $row=join("\t", @MetaT_id);
print OUT "Head\t$row\n";
foreach my $tmp1 (@Fun_id)
{
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (@MetaT_id)
        {
                if (exists $Fun2MetaT{$tmp1}{$tmp2})
                {
                        push @tmp, $Fun2MetaT{$tmp1}{$tmp2};
                }
                else
                {
                        push @tmp,"0"
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;



