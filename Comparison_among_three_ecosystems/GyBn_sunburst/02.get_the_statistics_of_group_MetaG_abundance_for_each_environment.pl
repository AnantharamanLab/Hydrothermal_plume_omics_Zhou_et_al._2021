#!/usr/bin/perl

use strict;
use warnings;

my %MetaG = (); #bin => row_mean => MetaG 10 times value
my @Head = ();my @Row_mean = ();
open IN, "MAG_average_coverage.row_mean.simple.txt";
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
		push @Row_mean, $key;
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

my %Group2Row_Mean_MetaG = (); # Group => Row_Mean => average MetaG 10 times value
foreach my $grp (sort keys %Microbial_group){
	foreach my $row_mean (@Row_mean){
		my $mean_metaG_value = 0; 
		foreach my $bin (sort keys %MetaG){
			if ($MAG_map{$bin}[0] eq  $grp){
				
				$mean_metaG_value += $MetaG{$bin}{$row_mean};
			}
		}
		
			
			$Group2Row_Mean_MetaG{$grp}{$row_mean} = $mean_metaG_value;
		
	}
}

#print out table
open OUT, ">MAG_average_coverage.Group2Row_Mean_MetaG.xls";
my $row=join("\t", @Row_mean);
print OUT "Head\t$row\n";
foreach my $tmp1 (sort keys %Group2Row_Mean_MetaG)
{
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (@Row_mean)
        {
                if (exists $Group2Row_Mean_MetaG{$tmp1}{$tmp2})
                {
                        push @tmp, $Group2Row_Mean_MetaG{$tmp1}{$tmp2};
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
open OUT, ">MAG_average_coverage.Group2Row_Mean_MetaG_each_bin.xls";
print OUT "group\tbin\tB.rowMean\tP.rowMean\tVF.rowMean\n";
foreach my $grp (sort keys %Microbial_group){	
	print OUT "$grp\t";
	foreach my $bin (sort keys %MetaG){		
		if ($MAG_map{$bin}[0] eq $grp){ 
			print OUT "$bin\t";
			my $row_mean_values = ""; my @Row_mean_values = ();
			foreach my $row_mean (@Row_mean){
				push @Row_mean_values, $MetaG{$bin}{$row_mean};
			}
			$row_mean_values = join("\t",@Row_mean_values);
			print OUT $row_mean_values."\n";
		}
	}
}
close OUT;
