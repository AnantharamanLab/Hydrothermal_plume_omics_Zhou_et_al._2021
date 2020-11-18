#!/usr/bin/perl

use strict;
use warnings;

my %MetaG = (); #bin => metaG_id => MetaG 10 times value
my @Head = ();my @MetaG_id = ();
open IN, "MetaG.10times.3_ecosystem.txt";
while (<IN>){
	chomp;
	if (/^Head/){
		my @tmp = split (/\t/,$_); @Head =@tmp;
	}else{
		my @tmp = split (/\t/,$_);
		my $bin = $tmp[0];
		for(my $i=1; $i<=$#tmp; $i++){
			$MetaG{$bin}{$Head[$i]} = $tmp[$i];
		}
	}
}
close IN;

foreach my $key (@Head){
	if ($key !~ "Head"){
		push @MetaG_id, $key;
	}
}

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

my %Group2MetaG_id_MetaG = (); # Group => MetaG_id => average MetaG 10 times value
foreach my $grp (sort keys %Microbial_group){
	foreach my $metaG_id (@MetaG_id){
		my $mean_metaG_value = 0; 
		foreach my $bin (sort keys %MetaG){
			if ($MAG_map{$bin}[0] eq  $grp){
				
				$mean_metaG_value += $MetaG{$bin}{$metaG_id};
			}
		}
		
			
			$Group2MetaG_id_MetaG{$grp}{$metaG_id} = $mean_metaG_value;
		
	}
}

#print out table
open OUT, ">MetaG.10times.3_ecosystem.Group2MetaG_id_MetaG.xls";
my $row=join("\t", @MetaG_id);
print OUT "Head\t$row\n";
foreach my $tmp1 (sort keys %Group2MetaG_id_MetaG)
{
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (@MetaG_id)
        {
                if (exists $Group2MetaG_id_MetaG{$tmp1}{$tmp2})
                {
                        push @tmp, $Group2MetaG_id_MetaG{$tmp1}{$tmp2};
                }
                else
                {
                        push @tmp,"0"
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;

#print out table2
open OUT, ">MetaG.10times.3_ecosystem.Group2MetaG_id_MetaG_each_bin.xls";
print OUT "group\tbin\tB.rowMean\tP.rowMean\tVF.rowMean\n";
foreach my $grp (sort keys %Microbial_group){	
	print OUT "$grp\t";
	foreach my $bin (sort keys %MetaG){		
		if ($MAG_map{$bin}[0] eq $grp){ 
			print OUT "$bin\t";
			my $metaG_id_values = ""; my @MetaG_id_values = ();
			foreach my $metaG_id (@MetaG_id){
				push @MetaG_id_values, $MetaG{$bin}{$metaG_id};
			}
			$metaG_id_values = join("\t",@MetaG_id_values);
			print OUT $metaG_id_values."\n";
		}
	}
}
close OUT;
