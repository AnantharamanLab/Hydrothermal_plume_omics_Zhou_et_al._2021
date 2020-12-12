#!/usr/bin/perl

use strict;
use warnings;

my %Fun2MetaT_abundance = (); # bin => metaT => 10times coverage

my @head = ();
open IN, "Fun2MetaT_abundance.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/,$_);
	if ($tmp[0] eq "Head"){
		@head = @tmp;
	}else{
		my ($Fun) = $tmp[0]; 
		for(my $i=1; $i<=$#head; $i++){
			$Fun2MetaT_abundance{$Fun}{$head[$i]} = $tmp[$i];
		}
	}	
}
close IN;

my %MetaT_id = (); my @MetaT_id = (); 
foreach my $key (@head){
	if ($key !~ /Head/){
		$MetaT_id{$key} = 1;
		push @MetaT_id, $key;
	}	
}

my %CymS_Fun_MetaT_HydrothermalPlume_coverage_Sig_result = (); #
my @CymS_Fun_MetaT_HydrothermalPlume_coverage_Sig_result = (); #
open IN, "CymS.Fun.MetaT.HydrothermalPlume.coverage.Sig.result.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/,$_);
	my $Fun = $tmp[0]; 
	my $output = $tmp[1]." | ".$tmp[0];
	$CymS_Fun_MetaT_HydrothermalPlume_coverage_Sig_result{$Fun} = $output;
	push @CymS_Fun_MetaT_HydrothermalPlume_coverage_Sig_result, $Fun;
}
close IN;

my @target_metaT_id_4_CymS_Fun_MetaT_HydrothermalPlume_coverage_Sig_result = qw /CymS.C.B	CymS.C.NBB.1	CymS.C.NBB.2	CymS.C.LRP.1	CymS.C.LRP.2	CymS.C.HRP.1	CymS.C.HRP.2/;


open OUT, ">CymS.Fun.MetaT.HydrothermalPlume.coverage.Sig.result.heatmap.table.txt";
my $row=join("\t", @target_metaT_id_4_CymS_Fun_MetaT_HydrothermalPlume_coverage_Sig_result);
print OUT "Head\t$row\n";
foreach my $tmp1 (@CymS_Fun_MetaT_HydrothermalPlume_coverage_Sig_result)
{
		print OUT $CymS_Fun_MetaT_HydrothermalPlume_coverage_Sig_result{$tmp1}."\t";
		my @tmp = ();
		foreach my $tmp2 (@target_metaT_id_4_CymS_Fun_MetaT_HydrothermalPlume_coverage_Sig_result)
		{
				if (exists $Fun2MetaT_abundance{$tmp1}{$tmp2})
				{
						push @tmp, $Fun2MetaT_abundance{$tmp1}{$tmp2};
				}
				else
				{
						push @tmp,"0"
				}
		}
		print OUT join("\t",@tmp)."\n";
}
close OUT;
