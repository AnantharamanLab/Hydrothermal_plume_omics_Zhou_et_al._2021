for i in `ls Cayman_gene_cat.gene*mapped.sorted.bam`
do
pileup.sh in=$i out=$i.pileup.out
done
