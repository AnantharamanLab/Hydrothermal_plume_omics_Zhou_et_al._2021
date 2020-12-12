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

#hmm result
my $hmm_result = "/mnt/storage7/zhouzhichao/BinProject/hydrothermal_plume_omics/MAG_analysis/do_hmm_annotation/All_faa.total.hmmscan.result.txt";
#custom_blastp result
my $blastp_result = "/mnt/storage7/zhouzhichao/BinProject/hydrothermal_plume_omics/MAG_analysis/custom_blastp/All_faa.total.custom_blastp.result.txt";
#KEGG_function result
my $KEGG_result = "/mnt/storage7/zhouzhichao/BinProject/hydrothermal_plume_omics/MAG_analysis/parse_KEGG_result/All_faa.total.KEGG_function.result.txt";
#dbCAN result
my $dbCAN_result = "/mnt/storage7/zhouzhichao/BinProject/hydrothermal_plume_omics/MAG_analysis/dbCAN/All_faa.total.dbCAN.result.txt";
#MEROPS result
my $MEROPS_result = "/mnt/storage7/zhouzhichao/BinProject/hydrothermal_plume_omics/MAG_analysis/MEROPS/Extracellular_analysis/All_faa.total.MEROPS.result.txt";


my %MAG_cov = (); # MetaG => Bin => coverage  100M reads normalized
my @head = ();
open IN, "/mnt/storage7/zhouzhichao/BinProject/hydrothermal_plume_omics/Cayman_Mapping.MetaT_new/MAG_average_coverage.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/,$_);
	if ($tmp[0] eq "Head"){
		@head = @tmp;
	}else{
		my $MetaG = $tmp[0] ; 
		for(my $i=1; $i<=$#head; $i++){
			$MAG_cov{$MetaG}{$head[$i]} = $tmp[$i];
		}
	}	
}
close IN;

my %MetaG_id = (); my @MetaG_id = (); 
open IN,"MetaG_list.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/,$_);
	if ($tmp[0] =~ /^Cym/){
		$MetaG_id{$tmp[0]} = 1;
		push @MetaG_id, $tmp[0];	
	}
}
close IN;

#-----------------------------------
my %Bin_id = (); # Store the Bin (aka. genome) names
my %Fun_result = (); # genome => function => number of hits
my @head_hmm_result = ();
open IN, $hmm_result;
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	if ($tmp[0] !~ /SZU/){
		@head_hmm_result = @tmp;
	}else{
		my ($genome) = $tmp[0] =~ /\_(SZU.+?)$/; #print "$genome\n";
		$Bin_id{$genome} = 1;
		for(my $i=1; $i<=$#head_hmm_result; $i++){
			$Fun_result{$genome}{$head_hmm_result[$i]} = $tmp[$i];
		}			
	}
}
close IN;

@head_hmm_result= ();
open IN, "/mnt/storage7/zhouzhichao/BinProject/hydrothermal_plume_omics/MAG_analysis/do_hmm_annotation/hmm_order.list";
while (<IN>){
	chomp;
	if (!/^#/){
		my @tmp = split (/\t/);
		my ($tmp2) = $tmp[2] =~ /^(.+?)\./;
		push @head_hmm_result, $tmp2;
	}
}
close IN;

my @head_blastp_result = qw /assA masD lactate_dehydrogenase phosphate_acetyltransferase acetate_kinase acetyl-CoA_synthetase aldehyde_dehydrogenase alcohol_dehydrogenase/;
	
open IN, $blastp_result;
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	if ($tmp[0] !~ /SZU/){
		@head_blastp_result = @tmp;
	}else{
		my ($genome) = $tmp[0] =~ /\_(SZU.+?)$/;
		for(my $i=1; $i<=$#head_blastp_result; $i++){
			$Fun_result{$genome}{$head_blastp_result[$i]} = $tmp[$i];
		}			
	}
}
close IN;

my @head_KEGG_result = ();
open IN, $KEGG_result;
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	if ($tmp[0] !~ /SZU/){
		@head_KEGG_result = @tmp;
	}else{
		my ($genome) = $tmp[0] =~ /\_(SZU.+?)$/;
		for(my $i=1; $i<=$#head_KEGG_result; $i++){
			$Fun_result{$genome}{$head_KEGG_result[$i]} = $tmp[$i];
		}			
	}
}
close IN;

my @head_dbCAN_result = ();
open IN, $dbCAN_result;
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	if ($tmp[0] !~ /SZU/){
		@head_dbCAN_result = @tmp;
	}else{
		my ($genome) = $tmp[0] =~ /\_(SZU.+?)$/;
		for(my $i=1; $i<=$#head_dbCAN_result; $i++){
			$Fun_result{$genome}{$head_dbCAN_result[$i]} = $tmp[$i];
		}			
	}
}
close IN;

my @head_MEROPS_result = ();
open IN, $MEROPS_result;
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	if ($tmp[0] !~ /SZU/){
		@head_MEROPS_result = @tmp;
	}else{
		my ($genome) = $tmp[0] =~ /\_(SZU.+?)$/;
		for(my $i=1; $i<=$#head_MEROPS_result; $i++){
			$Fun_result{$genome}{$head_MEROPS_result[$i]} = $tmp[$i];
		}			
	}
}
close IN;
#-----------------------------------

#-----------------------------------
my %Fun_id = (); # Store the function ids
my @Fun_id = ();
foreach my $key (@head_hmm_result){
	if ($key !~ /Bin|Head/){
		if (!exists $Fun_id{$key}){
			$Fun_id{$key} = 1;
			push @Fun_id,$key;
		}
	}
}

foreach my $key (@head_blastp_result){
	if ($key !~ /Bin|Head/){
		if (!exists $Fun_id{$key}){
			$Fun_id{$key} = 1;
			push @Fun_id,$key;
		}
	}
}

foreach my $key (@head_KEGG_result){
	if ($key !~ /Bin|Head/){
		if (!exists $Fun_id{$key}){
			$Fun_id{$key} = 1;
			push @Fun_id,$key;
		}
	}
}

foreach my $key (@head_dbCAN_result){
	if ($key !~ /Bin|Head/){
		if (!exists $Fun_id{$key}){
			$Fun_id{$key} = 1;
			push @Fun_id,$key;
		}
	}
}

foreach my $key (@head_MEROPS_result){
	if ($key !~ /Bin|Head/){
		if (!exists $Fun_id{$key}){
			$Fun_id{$key} = 1;
			push @Fun_id,$key;
		}
	}
}

#-----------------------------------

my %Fun2MetaG = (); #Fun => MetaG => Fun normalized abundance   #Bin is equal to genome in this script
foreach my $fun (sort keys %Fun_id){
	foreach my $MetaG (sort keys %MetaG_id){
		foreach my $Bin (sort keys %Bin_id){
			if ($Fun_result{$Bin}{$fun}){
				if ($MAG_cov{$MetaG}{$Bin}){
					$Fun2MetaG{$fun}{$MetaG} += $MAG_cov{$MetaG}{$Bin} * $Fun_result{$Bin}{$fun}; # %Fun_result: genome => function => number of hits
				}else{
					$Fun2MetaG{$fun}{$MetaG} += 0;
				}
			}
		}
	}	
}


open OUT, ">Fun2MetaG_abundance.txt";
my $row=join("\t", @MetaG_id);
print OUT "Head\t$row\n";
foreach my $tmp1 (@Fun_id)
{
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (@MetaG_id)
        {
                if (exists $Fun2MetaG{$tmp1}{$tmp2})
                {
                        push @tmp, $Fun2MetaG{$tmp1}{$tmp2};
                }
                else
                {
                        push @tmp,"0"
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;



