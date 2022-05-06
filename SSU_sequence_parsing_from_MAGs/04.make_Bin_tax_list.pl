#!/usr/bin/perl  -w

open IN,"<",  "checkm_result.txt" or die $!;
while (my $line = <IN>)
{
	if($line =~ /\_/)
	{
	my @tmp = split(/\s/,$line);
	my $fa ="";
	my $checkm = "";
	foreach (@tmp)
	{
		if ($_ =~ /scaf2bin/)
		{
			 $fa = $_;
		}
		if ($_ =~ /\_\_/ or $_ =~/root/)
		{
			 $checkm = $_;
		}
	}
	$hash{$fa}[0] = $line;		
	$hash{$fa}[1]= $checkm;
	}
}
close IN;

open IN, "<" , "checkm2tax";
while (my $line =<IN>)
{
	my @tmp = split(/\t/,$line);
	chomp $tmp[1];
	$hash2{$tmp[0]} = $tmp[1];
#	print $tmp[0]."\n";
}	
close IN;

foreach my $key (keys %hash)
{		
	if (exists $hash2{$hash{$key}[1]}){
		$hash{$key}[2] = $hash2{$hash{$key}[1]};
	}
}

open OUT, ">", "checkm_result_Bacteria.txt";
open OUT2, ">Bin_tax.list";
open OUT3, ">", "checkm_result_Archaea.txt";
foreach my $key (keys %hash)
{
	if ($hash{$key}[2] eq "Bacteria"){
		print OUT $hash{$key}[0];
		print OUT2 "$key\tBacteria\n";
		
	}elsif ($hash{$key}[2] eq "Archaea"){
		print OUT3 $hash{$key}[0];
		print OUT2 "$key\tArchaea\n";
	}else {
		print OUT2 "$key\troot\n";
	}

}
close OUT;
close OUT2;
close OUT3;
