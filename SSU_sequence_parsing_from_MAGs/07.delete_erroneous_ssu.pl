#!/usr/bin/perl

use strict;
use warnings;

my %Err16S = (); # scf_id => 1
open IN, "MAG_tax_16S_tax_comparing_table.txt";
while (<IN>){
	chomp;
	if (!/^Scaffold/){
		my @tmp = split (/\t/);
		my $scf_id = $tmp[0];
		if (defined $tmp[13] and $tmp[13] eq '0'){ # if the scaffold contains the erroness 16S fragment
			my $scf_id_new = $scf_id;
			if ($scf_id =~ /\#/){
				($scf_id_new) = $scf_id =~ /^(.+?)\-\#/;
			}
			$Err16S{$scf_id_new} = 1;
		}
	}
}
close IN;

`mkdir genomes/new_genomes`;
open IN, "ls genomes/*.fasta|";
while (<IN>){
	chomp;
	my $file = $_;
	my ($gn) = $file =~ /genomes\/(.+?)\.fasta/;
	my %Seq = _store_seq("$file");
	foreach my $head (sort keys %Seq){
		my $scf_id = $head;
		$scf_id =~ s/^>//g;
		$scf_id = $gn."\&\&".$scf_id;
		if (exists $Err16S{$scf_id}){
			delete $Seq{$head};
			#print "$head\n";
		}
	}
	
	open OUT, ">genomes/new_genomes/$gn.fasta";
	foreach my $head (sort keys %Seq){
		print OUT "$head\n$Seq{$head}\n";
	}
	close OUT;
}
close IN;

sub _store_seq{
	my $file = $_[0];
	my %Seq = (); my $head = "";
	open _IN, "$file";
	while (<_IN>){
		chomp;
		if (/>/){
			if (/\s/){
				($head) = $_ =~ /^(>.+?)\s/;
				$Seq{$head} = "";
			}else{
				($head) = $_ =~ /^(>.+?)$/;
				$Seq{$head} = "";
			}
		}else{
			$Seq{$head} .= $_;
		}
	}
	close _IN;
	return %Seq;
}



