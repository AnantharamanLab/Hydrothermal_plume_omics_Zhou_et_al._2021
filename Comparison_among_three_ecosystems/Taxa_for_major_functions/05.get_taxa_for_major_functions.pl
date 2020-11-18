#!/usr/bin/perl

use strict;
use warnings;

my %MAG_map = (); # bin => [0] Microbial group; [1] Genome Phylogeny; [2] Strain name
my %Microbial_group = ();
open IN,"MAG_info.txt";
while (<IN>){
	chomp;
	if (!/^#/){
		my @tmp = split (/\t/,$_);
		$MAG_map{$tmp[0]}[0] = $tmp[2]; $Microbial_group{$tmp[2]} = 1;
		$MAG_map{$tmp[0]}[1] = $tmp[4];
		$MAG_map{$tmp[0]}[2] = $tmp[5];
	}
}
close IN;

my %Fun_result= (); # fun_id => bin => no of hits
my @Fun_id = (); 
my @head = (); #Store the head line
open IN, "Fun_result.txt";
while (<IN>){
	chomp;
	if (/^Head/){
		my @tmp = split (/\t/); @head = @tmp;
		foreach my $key (@tmp){
			if ($key !~ /Head/){
				push @Fun_id, $key;
			}
		}		
	}else{
		my @tmp = split (/\t/);
		for(my $i=1; $i<=$#tmp; $i++){
			$Fun_result{$head[$i]}{$tmp[0]} = $tmp[$i];
		}
	}
}
close IN;

my %MetaG2Cov= (); # metaG => bin => 10 times coverage value
my @MetaG_id = (); 
my @head2 = (); #Store the head line
open IN, "MetaG.10times_coverage.txt";
while (<IN>){
	chomp;
	if (/^Head/){
		my @tmp = split (/\t/); @head2 = @tmp;
		foreach my $key (@tmp){
			if ($key !~ /Head/){
				push @MetaG_id, $key;
			}
		}		
	}else{
		my @tmp = split (/\t/);
		for(my $i=1; $i<=$#tmp; $i++){
			$MetaG2Cov{$head2[$i]}{$tmp[0]} = $tmp[$i];
		}
	}
}
close IN;

my %Major_fun = (); #B.Lau.mean => TIGR00014
my %Fun_id_this_time = (); my @Fun_id_this_time = ();
open IN, "Major_functions.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	$Major_fun{$tmp[0]} = $tmp[1];
	$Fun_id_this_time{$tmp[1]} = 1;
	push @Fun_id_this_time,$tmp[1];
}
close IN;


my %Major_fun_micro_grp = (); #B.Lau.mean => fun_id => microbial group => percentage;
foreach my $metaG (sort keys %Major_fun){
	foreach my $fun_id (sort keys %Fun_id_this_time){
		my $total_fun_id_cov = 0; # the coverage value of fun trait for certain function id; coverage value = fun trait hits in the bin multiplied by bin coverage 
		my %Micro_grp_fun_id_cov = (); #the coverage value of fun trait in certain microbial groups
		
			foreach my $bin (sort keys %MAG_map){
				if ($Fun_result{$fun_id}{$bin}){
					if ($MetaG2Cov{$metaG}{$bin}){
						$total_fun_id_cov += $MetaG2Cov{$metaG}{$bin} * $Fun_result{$fun_id}{$bin};
					}else{
						$total_fun_id_cov += 0;
					}						
				}
				foreach my $micro_grp (sort keys %Microbial_group){
					if ($micro_grp eq $MAG_map{$bin}[0]){
						if ($Fun_result{$fun_id}{$bin}){
							if ($MetaG2Cov{$metaG}{$bin}){
								$Micro_grp_fun_id_cov{$micro_grp} += $MetaG2Cov{$metaG}{$bin} * $Fun_result{$fun_id}{$bin};
							}else{
								$Micro_grp_fun_id_cov{$micro_grp} += 0;
							}						
						}
					}
				}
			}
		
		foreach my $micro_grp (sort keys %Micro_grp_fun_id_cov){
			if ($Micro_grp_fun_id_cov{$micro_grp}){
				$Major_fun_micro_grp{$metaG}{$fun_id}{$micro_grp} = $Micro_grp_fun_id_cov{$micro_grp} / $total_fun_id_cov;
			}else{
				$Major_fun_micro_grp{$metaG}{$fun_id}{$micro_grp} = 0;
			}
		}		
	}	
}

foreach my $metaG (sort keys %Major_fun){
	open OUT, ">Major_fun_micro_grp.$metaG.txt";
	my $row=join("\t", @Fun_id_this_time);
	print OUT "Head\t$row\n";
	foreach my $tmp1 (sort keys %Microbial_group)
	{
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (@Fun_id_this_time)
        {
                if (exists $Major_fun_micro_grp{$metaG}{$tmp2}{$tmp1})
                {
                        push @tmp, $Major_fun_micro_grp{$metaG}{$tmp2}{$tmp1};
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




