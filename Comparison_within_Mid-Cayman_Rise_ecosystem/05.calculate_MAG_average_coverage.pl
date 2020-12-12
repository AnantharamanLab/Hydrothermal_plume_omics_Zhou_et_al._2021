#!/usr/bin/perl

use strict;
use warnings;
use Statistics::Descriptive;

my %map = (); # DNA and cDNA reads maps
open IN, "Transcriptom_map.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	$map{$tmp[3]}[0] = $tmp[0];   # SRR [0] => name
	$map{$tmp[3]}[1] = $tmp[1];   # SRR [1] => cDNA or DNA
	$map{$tmp[3]}[2] = $tmp[2];   # SRR [2] => PAIRED or SINGLE
	$map{$tmp[3]}[3] = $tmp[4];   # SRR [3] => reads numbers
}
close IN;

my %h = (); my @h_head = ();  my @h_head_num = (); my %Bin = ();
open IN, "Cayman_genome_cat.genome.depth.txt";
while (<IN>){
	chomp;
	if (/^contigName/){
		my @tmp = split (/\t/);@h_head = @tmp;		
		for(my $i=0; $i<=$#h_head; $i++){
			if ($h_head[$i] =~ /.sorted.bam$/){
					push @h_head_num, $i;
			}			
		}
	}else{
		my @tmp = split (/\t/);
		my ($bin) = $tmp[0] =~ /^(.+?)\_/;
		$Bin{$bin} = 1;
		foreach my $i (@h_head_num){
			if (!exists $h{$h_head[$i]}{$bin}){
				$h{$h_head[$i]}{$bin} = $tmp[$i];
			}else{
				$h{$h_head[$i]}{$bin} .= "\t".$tmp[$i];
			}
		}
	}
}
close IN;

foreach my $i (@h_head_num){
	foreach my $bin (sort keys %Bin){
		my @tmp = split (/\t/, $h{$h_head[$i]}{$bin});
		my $stat = Statistics::Descriptive::Full->new();
		$stat->add_data(\@tmp);
		my $mean = $stat->mean();
		my ($name) = $h_head[$i] =~ /\_cat\.genome\_(.+?)\_mapped\.sorted\.bam/; 
		my $reads_num = 0;
		foreach my $key (sort keys %map){
			if ($map{$key}[0] eq $name){
				$reads_num = $map{$key}[3];
			}
		}
		$h{$h_head[$i]}{$bin} = $mean * (100000000 / $reads_num ); #100M reads
	}	
}

#print table
open OUT, ">MAG_average_coverage.txt";
my $row=join("\t", sort keys %Bin);
print OUT "Head\t$row\n";
foreach my $tmp1 (sort keys %h)
{
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (sort keys %Bin)
        {
                if (exists $h{$tmp1}{$tmp2})
                {
                        push @tmp, $h{$tmp1}{$tmp2};
                }
                else
                {
                        push @tmp,"0"
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;

open OUT, ">tmp.MAG_average_coverage.txt";
open IN, "MAG_average_coverage.txt";
while (<IN>){
	chomp;
	if (!/sorted\.bam/){
		print OUT "$_\n";
	}else{
		my $tmp_out = $_;
		$tmp_out =~ s/Cayman\_genome\_cat\.genome_//g;
		$tmp_out =~ s/\_mapped\.sorted\.bam//g;
		print OUT "$tmp_out\n";
	}
}
close IN;
close OUT;

`rm MAG_average_coverage.txt`;
`mv tmp.MAG_average_coverage.txt MAG_average_coverage.txt`;
