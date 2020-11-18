for i in `ls *all.gene*mapped.sorted.bam`
do
pileup.sh in=$i out=$i.pileup.out
done
