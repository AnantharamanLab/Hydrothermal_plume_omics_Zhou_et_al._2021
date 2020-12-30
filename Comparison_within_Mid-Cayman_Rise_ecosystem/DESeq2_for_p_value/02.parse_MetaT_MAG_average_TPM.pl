#!/usr/bin/perl

use strict;
use warnings;
	
my %MetaT_TPM = (); # gene => MetaT => TPM * 100  
my @head = ();
my %gene_id = ();
my %Bin_id = ();
open IN, "MetaT.TPM.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/,$_);
	if ($tmp[0] eq "Head"){
		@head = @tmp;
	}else{
		my $gene = $tmp[0]; $gene_id{$gene} = 1; 
		my ($bin) = $gene =~ /^.+?\_(SZU.+?)\|/; $Bin_id{$bin} = 1; 
		for(my $i=1; $i<=$#head; $i++){
			$MetaT_TPM{$gene}{$head[$i]} = $tmp[$i] * 100;
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

my %Bin2MetaT = (); # bin => metaT => TPM
foreach my $gene (sort keys %gene_id){
	foreach my $MetaT (sort keys %MetaT_id){
		my ($bin) = $gene =~ /^.+?\_(SZU.+?)\|/; 
		if (!exists $Bin2MetaT{$bin}{$MetaT}){
			$Bin2MetaT{$bin}{$MetaT} = $MetaT_TPM{$gene}{$MetaT};
		}else{
			$Bin2MetaT{$bin}{$MetaT} .= "\t".$MetaT_TPM{$gene}{$MetaT};
		}
	}
}
	
my %Bin2MetaT2 = (); #The new %Bin2MetaT2 that store the average value of all TPM
foreach my $bin (sort keys %Bin2MetaT){
	foreach my $MetaT (sort keys %MetaT_id){
		my @tmp = split (/\t/,$Bin2MetaT{$bin}{$MetaT});
		my $tmp_sum = 0; 
		foreach my $item (@tmp){
			$tmp_sum += $item;
		}
		$Bin2MetaT2{$bin}{$MetaT} = $tmp_sum / (scalar @tmp);
	}	
}

open OUT, ">Bin2MetaT_abundance.txt";
my $row=join("\t", @MetaT_id);
print OUT "Head\t$row\n";
foreach my $tmp1 (sort keys %Bin_id)
{
		print OUT $tmp1."\t";
		my @tmp = ();
		foreach my $tmp2 (@MetaT_id)
		{
				if (exists $Bin2MetaT2{$tmp1}{$tmp2})
				{
						push @tmp, $Bin2MetaT2{$tmp1}{$tmp2};
				}
				else
				{
						push @tmp,"0"
				}
		}
		print OUT join("\t",@tmp)."\n";
}
close OUT;








