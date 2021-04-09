/slowdata/archive/anaconda2/bin/mmseqs easy-cluster All_plume_phage_drep95_rep_seq.fasta All_plume_phage_drep95_looseCluster2 ./temp --min-seq-id 0.25 -c 0.5 -s 7.5 -e 0.001 --threads 10

cat All_plume_phage_drep95_looseCluster2_all_seqs.fasta | sed 's/Axial.Plume~~/Axial~~/g' | sed 's/Axial.Seawater~~/Axial~~/g' > All_plume_phage_drep95_looseCluster2_all_seqs_mdf.fasta


# previously we use --min-seq-id 0.5 -c 0.6 -s 7.5 -e 0.001

#--min-seq-id (sequence identity cutoff)

#-c (coverage cutoff)

#-s (sensitivity, usually just stay at the max of 7.5)

#-e (evalue)