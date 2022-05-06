#multi-core version
N=150
for i in `ls genomes/*.fasta`
do
	k=$(basename $i .fasta)
        ((j=j%N)); ((j++==0)) && wait
        checkm ssu_finder -x .fasta  genomes/$k.fasta  ./genomes ./checkm_result/ssu_finder/${k}_result -t 40 &
done

