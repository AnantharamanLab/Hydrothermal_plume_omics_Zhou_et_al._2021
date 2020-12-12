#!/usr/bin/perl

use strict;
use warnings;

my %CymD_100times_MAG_MetaT_coverage = (); # bin => metaT => 100 times coverage

my @head = ();
open IN, "CymD.100times.MAG.MetaT.coverage.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/,$_);
	if ($tmp[0] eq "Head"){
		@head = @tmp;
	}else{
		my ($bin) = $tmp[0]; 
		for(my $i=1; $i<=$#head; $i++){
			$CymD_100times_MAG_MetaT_coverage{$bin}{$head[$i]} = $tmp[$i];
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

my %CymD_MAG_MetaT_HydrothermalPlume_coverage_Sig_result = (); #
my @CymD_MAG_MetaT_HydrothermalPlume_coverage_Sig_result = (); #
open IN, "CymD.MAG.MetaT.HydrothermalPlume.coverage.Sig.result.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/,$_);
	my $bin = $tmp[0]; 
	my $tax = $tmp[1]." | ".$tmp[3]." | ".$tmp[5];
	$CymD_MAG_MetaT_HydrothermalPlume_coverage_Sig_result{$bin} = $tax;
	push @CymD_MAG_MetaT_HydrothermalPlume_coverage_Sig_result, $bin;
}
close IN;

my @target_metaT_id_4_CymD_MAG_MetaT_HydrothermalPlume_coverage_Sig_result = qw /CymD.C.LB	CymD.C.HB	CymD.C.LRP	CymD.C.MRP	CymD.C.HRP/;


open OUT, ">CymD.MAG.MetaT.HydrothermalPlume.coverage.Sig.result.heatmap.table.txt";
my $row=join("\t", @target_metaT_id_4_CymD_MAG_MetaT_HydrothermalPlume_coverage_Sig_result);
print OUT "Head\t$row\n";
foreach my $tmp1 (@CymD_MAG_MetaT_HydrothermalPlume_coverage_Sig_result)
{
		print OUT $CymD_MAG_MetaT_HydrothermalPlume_coverage_Sig_result{$tmp1}."\t";
		my @tmp = ();
		foreach my $tmp2 (@target_metaT_id_4_CymD_MAG_MetaT_HydrothermalPlume_coverage_Sig_result)
		{
				if (exists $CymD_100times_MAG_MetaT_coverage{$tmp1}{$tmp2})
				{
						push @tmp, $CymD_100times_MAG_MetaT_coverage{$tmp1}{$tmp2};
				}
				else
				{
						push @tmp,"0"
				}
		}
		print OUT join("\t",@tmp)."\n";
}
close OUT;
