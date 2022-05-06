#!/usr/bin/perl  -w

my %hash = ();
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
                if ($_ =~ /\_\_/ or $_ =~/root/)
                {
                         $checkm = $_;
                }
        }
        $hash{$checkm} = "";
        }
}
close IN;

open OUT, ">checkm2tax";
foreach my $key (sort keys %hash){
	print OUT "$key\n";
}
close OUT;
