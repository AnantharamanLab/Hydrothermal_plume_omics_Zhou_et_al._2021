#!/usr/bin/perl

use strict;
use warnings;

my %Fun2MetaG_abundance = (); # bin => metaG => 10times coverage

my @head = ();
open IN, "Fun2MetaG_abundance.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/,$_);
	if ($tmp[0] eq "Head"){
		@head = @tmp;
	}else{
		my ($Fun) = $tmp[0]; 
		for(my $i=1; $i<=$#head; $i++){
			$Fun2MetaG_abundance{$Fun}{$head[$i]} = $tmp[$i];
		}
	}	
}
close IN;

my %MetaG_id = (); my @MetaG_id = (); 
foreach my $key (@head){
	if ($key !~ /Head/){
		$MetaG_id{$key} = 1;
		push @MetaG_id, $key;
	}	
}

my %CymS_Fun_MetaG_HydrothermalPlume_coverage_Sig_result = (); #
my @CymS_Fun_MetaG_HydrothermalPlume_coverage_Sig_result = (); #
open IN, "CymS.Fun.MetaG.HydrothermalPlume.coverage.Sig.result.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/,$_);
	my $Fun = $tmp[0]; 
	my $output = $tmp[1]." | ".$tmp[0];
	$CymS_Fun_MetaG_HydrothermalPlume_coverage_Sig_result{$Fun} = $output;
	push @CymS_Fun_MetaG_HydrothermalPlume_coverage_Sig_result, $Fun;
}
close IN;

my @target_metaG_id_4_CymS_Fun_MetaG_HydrothermalPlume_coverage_Sig_result = qw /CymS.D.B CymS.D.LRP.1 CymS.D.HRP.1/;


open OUT, ">CymS.Fun.MetaG.HydrothermalPlume.coverage.Sig.result.heatmap.table.txt";
my $row=join("\t", @target_metaG_id_4_CymS_Fun_MetaG_HydrothermalPlume_coverage_Sig_result);
print OUT "Head\t$row\n";
foreach my $tmp1 (@CymS_Fun_MetaG_HydrothermalPlume_coverage_Sig_result)
{
		print OUT $CymS_Fun_MetaG_HydrothermalPlume_coverage_Sig_result{$tmp1}."\t";
		my @tmp = ();
		foreach my $tmp2 (@target_metaG_id_4_CymS_Fun_MetaG_HydrothermalPlume_coverage_Sig_result)
		{
				if (exists $Fun2MetaG_abundance{$tmp1}{$tmp2})
				{
						push @tmp, $Fun2MetaG_abundance{$tmp1}{$tmp2};
				}
				else
				{
						push @tmp,"0"
				}
		}
		print OUT join("\t",@tmp)."\n";
}
close OUT;
