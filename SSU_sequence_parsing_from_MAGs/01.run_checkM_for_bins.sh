#checkm lineage_wf -x .fasta ./genomes ./checkm_result -f checkm_result.txt -t 40

mkdir tmpdir
checkm lineage_wf -x .fasta -t 40 --pplacer_threads 40 --tmpdir ./tmpdir/ ./genomes  ./checkm_result -f ./checkm_result.txt
rm -r tmpdir
