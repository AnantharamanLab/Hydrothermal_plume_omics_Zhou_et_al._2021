#!/usr/bin/perl

use strict;
use warnings;

my %CymS_100times_MAG_MetaT_coverage = (); # bin => metaT => 10times coverage

my @head = ();
open IN, "CymS.100times.MAG.MetaT.coverage.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/,$_);
	if ($tmp[0] eq "Head"){
		@head = @tmp;
	}else{
		my ($bin) = $tmp[0]; 
		for(my $i=1; $i<=$#head; $i++){
			$CymS_100times_MAG_MetaT_coverage{$bin}{$head[$i]} = $tmp[$i];
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

my %CymS_MAG_MetaT_HydrothermalPlume_coverage_Sig_result = (); #
my @CymS_MAG_MetaT_HydrothermalPlume_coverage_Sig_result = (); #
open IN, "CymS.MAG.MetaT.HydrothermalPlume.coverage.Sig.result.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/,$_);
	my $bin = $tmp[0]; 
	my $tax = $tmp[1]." | ".$tmp[3]." | ".$tmp[5];
	$CymS_MAG_MetaT_HydrothermalPlume_coverage_Sig_result{$bin} = $tax;
	push @CymS_MAG_MetaT_HydrothermalPlume_coverage_Sig_result, $bin;
}
close IN;

my @target_metaT_id_4_CymS_MAG_MetaT_HydrothermalPlume_coverage_Sig_result = qw /CymS.C.B	CymS.C.NBB.1	CymS.C.NBB.2	CymS.C.LRP.1	CymS.C.LRP.2	CymS.C.HRP.1	CymS.C.HRP.2/;


open OUT, ">CymS.MAG.MetaT.HydrothermalPlume.coverage.Sig.result.heatmap.table.txt";
my $row=join("\t", @target_metaT_id_4_CymS_MAG_MetaT_HydrothermalPlume_coverage_Sig_result);
print OUT "Head\t$row\n";
foreach my $tmp1 (@CymS_MAG_MetaT_HydrothermalPlume_coverage_Sig_result)
{
		print OUT $CymS_MAG_MetaT_HydrothermalPlume_coverage_Sig_result{$tmp1}."\t";
		my @tmp = ();
		foreach my $tmp2 (@target_metaT_id_4_CymS_MAG_MetaT_HydrothermalPlume_coverage_Sig_result)
		{
				if (exists $CymS_100times_MAG_MetaT_coverage{$tmp1}{$tmp2})
				{
						push @tmp, $CymS_100times_MAG_MetaT_coverage{$tmp1}{$tmp2};
				}
				else
				{
						push @tmp,"0"
				}
		}
		print OUT join("\t",@tmp)."\n";
}
close OUT;
