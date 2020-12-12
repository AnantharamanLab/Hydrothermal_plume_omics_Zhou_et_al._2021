#!/usr/bin/perl

use strict;
use warnings;

my %CymS_10times_transpose_MAG_MetaG_coverage = (); # bin => metaG => 10times coverage

my @head = ();
open IN, "CymS.10times.transpose.MAG.MetaG.coverage.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/,$_);
	if ($tmp[0] eq "Head"){
		@head = @tmp;
	}else{
		my ($bin) = $tmp[0]; 
		for(my $i=1; $i<=$#head; $i++){
			$CymS_10times_transpose_MAG_MetaG_coverage{$bin}{$head[$i]} = $tmp[$i];
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

my %CymS_MAG_MetaG_HydrothermalPlume_coverage_Sig_result = (); #
my @CymS_MAG_MetaG_HydrothermalPlume_coverage_Sig_result = (); #
open IN, "CymS.MAG.MetaG.HydrothermalPlume.coverage.Sig.result.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/,$_);
	my $bin = $tmp[0]; 
	my $tax = $tmp[1]." | ".$tmp[3]." | ".$tmp[5];
	$CymS_MAG_MetaG_HydrothermalPlume_coverage_Sig_result{$bin} = $tax;
	push @CymS_MAG_MetaG_HydrothermalPlume_coverage_Sig_result, $bin;
}
close IN;

my @target_metaG_id_4_CymS_MAG_MetaG_HydrothermalPlume_coverage_Sig_result = qw /CymS.D.B	CymS.D.LRP.1 CymS.D.HRP.1/;


open OUT, ">CymS.MAG.MetaG.HydrothermalPlume.coverage.Sig.result.heatmap.table.txt";
my $row=join("\t", @target_metaG_id_4_CymS_MAG_MetaG_HydrothermalPlume_coverage_Sig_result);
print OUT "Head\t$row\n";
foreach my $tmp1 (@CymS_MAG_MetaG_HydrothermalPlume_coverage_Sig_result)
{
		print OUT $CymS_MAG_MetaG_HydrothermalPlume_coverage_Sig_result{$tmp1}."\t";
		my @tmp = ();
		foreach my $tmp2 (@target_metaG_id_4_CymS_MAG_MetaG_HydrothermalPlume_coverage_Sig_result)
		{
				if (exists $CymS_10times_transpose_MAG_MetaG_coverage{$tmp1}{$tmp2})
				{
						push @tmp, $CymS_10times_transpose_MAG_MetaG_coverage{$tmp1}{$tmp2};
				}
				else
				{
						push @tmp,"0"
				}
		}
		print OUT join("\t",@tmp)."\n";
}
close OUT;
