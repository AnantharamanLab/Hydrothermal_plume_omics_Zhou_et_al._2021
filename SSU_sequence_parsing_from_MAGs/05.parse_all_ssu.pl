#!/usr/bin/perl -w 

use strict;
use warnings;

my %tax = ();
open IN, "Bin_tax.list";
while (<IN>){
	chomp;
	my @tmp = split (/\t/,$_);
	$tax{$tmp[0]} = $tmp[1];
}
close IN;

my $head_scf = "Scaffold ID\tHMM\ti-Evalue\tStart hit\tEnd hit\t16S/18S gene length\tRev. Complement\tSequence length"; # store the head line of scf info
my %Total_scf = (); # $scf_id => whole line
my %out = ();
open IN, "ls checkm_result/ssu_finder/*/ssu_summary.tsv | ";
while (<IN>){
	chomp;
	my ($fa) = $_ =~ /ssu_finder\/(.+?)\_result/;
	my $dir = $_; $dir =~ s/\/ssu_summary\.tsv//g;
	my %hit = (); # $scf_id => whole line 
	open INN, "$_";
	OUTTER1: while (<INN>){
		chomp;
		if (!/^Bin/){
			my @tmp = split (/\t/,$_) or last OUTTER1;
			my $scf_id = $fa."\&\&".$tmp[1];
			$hit{$scf_id} = $scf_id."\t".$tmp[2]."\t".$tmp[3]."\t".$tmp[4]."\t".$tmp[5]."\t".$tmp[6]."\t".$tmp[7]."\t".$tmp[8];
		}
	}
	close INN;
	if (%hit){
		#Store ssu.fna
		my %ssu = ();my $head = "";
		open _IN, "$dir/ssu.fna" or die ; 
		while (<_IN>){
			chomp;
			if (/>/){
				$head = $_;
				$ssu{$head} = "";
			}else{
				$ssu{$head} .= $_;
			}
		}
		close _IN; 
		
		%out = (%out, %ssu);
		%Total_scf = (%Total_scf, %hit);
	}
}
close IN;

open OUT, ">Bin_ssu_all.fa";
foreach my $key (sort keys %out){
	print OUT "$key\n$out{$key}\n";
}
close OUT;

open OUT, ">Total_scf_info.txt";
print OUT  $head_scf."\n";
foreach my $scf_id (sort keys %Total_scf){
	print OUT "$Total_scf{$scf_id}\n";
}
close OUT;
