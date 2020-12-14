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
#hmm iron_and_magnesium_oxidation dir; result files are ended with .hmmscan.result.txt
my $hmm_iron_and_magnesium_oxidation_result = "/mnt/storage7/zhouzhichao/BinProject/hydrothermal_plume_omics/MAG_analysis/do_hmm_annotation_iron_magnesium/All_faa.total.hmmscan.result.txt";
#custom_blastp result
my $blastp_result = "/mnt/storage7/zhouzhichao/BinProject/hydrothermal_plume_omics/MAG_analysis/custom_blastp/All_faa.total.custom_blastp.result.txt";
#KEGG_function result
my $KEGG_result = "/mnt/storage7/zhouzhichao/BinProject/hydrothermal_plume_omics/MAG_analysis/parse_KEGG_result/All_faa.total.KEGG_function.result.txt";
#dbCAN result
my $dbCAN_result = "/mnt/storage7/zhouzhichao/BinProject/hydrothermal_plume_omics/MAG_analysis/dbCAN/All_faa.total.dbCAN.result.txt";
#MEROPS result
my $MEROPS_result = "/mnt/storage7/zhouzhichao/BinProject/hydrothermal_plume_omics/MAG_analysis/MEROPS/Extracellular_analysis/All_faa.total.MEROPS.result.txt";

open INN_, "ls /mnt/storage7/zhouzhichao/BinProject/hydrothermal_plume_omics/Sep_Mapping.MetaT/*.MAG_average_coverage.txt |";
while (<INN_>){
	chomp;
	my $file = $_;
	my ($env) = $file =~ /Sep_Mapping\.MetaT\/(.+?)\.MAG\_average\_coverage\.txt/;
	my %MAG_cov = (); # MetaG => Bin => coverage 100M reads normalized
	my @head = ();
	my %Bin_id = ();
	my %MetaG_id = ();
	my @MetaG_id = ();
	open IN, $file;
	while (<IN>){
		chomp;
		my @tmp = split (/\t/,$_);
		if ($tmp[0] eq "Head"){
			@head = @tmp;
		}else{
			my ($MetaG) = $tmp[0] =~ /genome\_(.+?)\_mapped/; $MetaG_id{$MetaG} = 1; push @MetaG_id, $MetaG;
			for(my $i=1; $i<=$#head; $i++){
				$MAG_cov{$MetaG}{$head[$i]} = $tmp[$i];
			}
		}	
	}
	close IN;
=pod
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
=cut
	my %FunMap = (); # Fun1 (target id) => Fun2 (electron donor or acceptor), for example: TIGR02066 => Electron donor|S0	
	open IN, "FunMap.list";
	while (<IN>){
		chomp;
		my @tmp = split (/\t/,$_);
		$FunMap{$tmp[1]} = $tmp[0];
	}
	close IN;


	#-----------------------------------
	#my %Bin_id = (); # Store the Bin (aka. genome) names
	my %Fun_result = (); # genome => function => number of hits
	my %Fun2_id = (); # Store the fun2 ids
	my @Fun2_id = ();
	my %Fun2_result = (); # genome => fun2 =>  number of hits
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
				my $fun = $head_hmm_result[$i];
				$Fun_result{$genome}{$fun} = $tmp[$i];
				if (exists $FunMap{$fun}){
					my $fun2 = $FunMap{$fun}; 
					if (!exists $Fun2_id{$fun2}){
						push @Fun2_id,$fun2;
					}
					$Fun2_id{$fun2} = 1; 
					$Fun2_result{$genome}{$fun2} += $tmp[$i];
				}
				
			}			
		}
	}
	close IN;
	
	open IN, $hmm_iron_and_magnesium_oxidation_result;
	while (<IN>){
		chomp;
		my @tmp = split (/\t/);
		if ($tmp[0] !~ /SZU/){
			@head_hmm_result = @tmp;
		}else{
			my ($genome) = $tmp[0] =~ /\_(SZU.+?)$/; #print "$genome\n";
			$Bin_id{$genome} = 1;
			for(my $i=1; $i<=$#head_hmm_result; $i++){
				my $fun = $head_hmm_result[$i];
				$Fun_result{$genome}{$fun} = $tmp[$i];
				if (exists $FunMap{$fun}){
					my $fun2 = $FunMap{$fun}; 
					if (!exists $Fun2_id{$fun2}){
						push @Fun2_id,$fun2;
					}
					$Fun2_id{$fun2} = 1; 
					$Fun2_result{$genome}{$fun2} += $tmp[$i];
				}
				
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
	
	open IN, "/mnt/storage7/zhouzhichao/BinProject/hydrothermal_plume_omics/MAG_analysis/do_hmm_annotation_iron_magnesium/hmm_order.list";
	while (<IN>){
		chomp;
		if (!/^#/){
			my @tmp = split (/\t/);
			my ($tmp2) = $tmp[2] =~ /^(.+?)\./;
			push @head_hmm_result, $tmp2;
		}
	}
	close IN;	

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


	#-----------------------------------

	my %Fun2_2_MetaG = (); #fun2_id => MetaG => Fun normalized abundance.  #Bin is equal to genome in this script. Fun2 (electron donor or acceptor) is the one with the ";".
	my %Fun2_2_MetaG2bin = (); #fun2_id => MetaG => bin => Fun normalized abundance.
	foreach my $fun2_id (sort keys %Fun2_id){
		foreach my $MetaG (sort keys %MetaG_id){
			foreach my $Bin (sort keys %Bin_id){
				my $fun_norm_abund  = 0; # for the result value of %Fun2_2_MetaG2bin; after each cycle, $fun_norm_abund for each bin will be reset to zero.
				foreach my $fun (sort keys %FunMap){
					if ($FunMap{$fun} eq $fun2_id){
						if ($Fun_result{$Bin}{$fun}){
							if ($MAG_cov{$MetaG}{$Bin}){
								$Fun2_2_MetaG{$fun2_id}{$MetaG} += $MAG_cov{$MetaG}{$Bin} * $Fun_result{$Bin}{$fun}; # %Fun_result: genome => function => number of hits
								$fun_norm_abund += $MAG_cov{$MetaG}{$Bin} * $Fun_result{$Bin}{$fun};
							}else{
								$Fun2_2_MetaG{$fun2_id}{$MetaG} += 0;
								$fun_norm_abund += 0;
							}
						}
					}				
				}
				$Fun2_2_MetaG2bin{$fun2_id}{$MetaG}{$Bin} = $fun_norm_abund;
			}
		}	
	}

	my %Fun2_2_MetaG_new = ();  #fun2_id(new) => MetaG => Fun normalized abundance. Some fun2_id have two ids which are connected by ";". 
	my %Fun2_id_new = (); #fun2_id (new) => 1;
	foreach my $fun2 (sort keys %Fun2_id){ 
		foreach my $MetaG (sort keys %MetaG_id){
			if ($fun2 =~ /\;/){
				my @tmp = split (/\;/,$fun2);
				my $fun2_1 = $tmp[0]; my $fun2_2 = $tmp[1]; $Fun2_id_new{$fun2_1} = 1; $Fun2_id_new{$fun2_2} = 1; 				
				if ($Fun2_2_MetaG{$fun2}{$MetaG}){
					$Fun2_2_MetaG_new{$fun2_1}{$MetaG} += $Fun2_2_MetaG{$fun2}{$MetaG};
					$Fun2_2_MetaG_new{$fun2_2}{$MetaG} += $Fun2_2_MetaG{$fun2}{$MetaG};
				}else{
					$Fun2_2_MetaG_new{$fun2_1}{$MetaG} += 0; $Fun2_2_MetaG_new{$fun2_2}{$MetaG} +=0;
				}
			}else{
				$Fun2_id_new{$fun2} = 1; 
				if ($Fun2_2_MetaG{$fun2}{$MetaG}){
					$Fun2_2_MetaG_new{$fun2}{$MetaG} += $Fun2_2_MetaG{$fun2}{$MetaG};
				}else{
					$Fun2_2_MetaG_new{$fun2}{$MetaG} += 0;
				}
			}		
		}	
	}

	my %Fun2_id_new2target_num = (); #Fun2 (new) => the number of target ids;
	foreach my $fun2_id_new (sort keys %Fun2_id_new){
		my $target_num = 0;
		foreach my $target_id (sort keys %FunMap){
			if ($FunMap{$target_id} =~ /$fun2_id_new/){
				$target_num++;
			}
		}
		$Fun2_id_new2target_num{$fun2_id_new} = $target_num;
	}

	#use the target num to divide the Fun normalized abundance of each fun2_id_new and each MetaG
	foreach my $fun2_id_new (sort keys %Fun2_id_new){ 
		foreach my $MetaG (sort keys %MetaG_id){
			if ($Fun2_id_new2target_num{$fun2_id_new} and $Fun2_2_MetaG_new{$fun2_id_new}{$MetaG}){
				$Fun2_2_MetaG_new{$fun2_id_new}{$MetaG} = $Fun2_2_MetaG_new{$fun2_id_new}{$MetaG} / $Fun2_id_new2target_num{$fun2_id_new};
			}
		}
	}

	open OUT, ">$env.Electron.Fun2MetaG_abundance.txt";
	my $row=join("\t", @MetaG_id);
	print OUT "Head\t$row\n";
	foreach my $tmp1 (sort keys %Fun2_id_new)
	{
			print OUT $tmp1."\t";
			my @tmp = ();
			foreach my $tmp2 (@MetaG_id)
			{
					if (exists $Fun2_2_MetaG_new{$tmp1}{$tmp2})
					{
							push @tmp, $Fun2_2_MetaG_new{$tmp1}{$tmp2};
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
close INN_;



