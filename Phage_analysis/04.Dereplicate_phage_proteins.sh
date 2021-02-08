/slowdata/archive/anaconda2/bin/mmseqs easy-cluster All_plume_phage_proteins.faa All_plume_phage_drep ./temp --min-seq-id 0.95 -c 0.8 -s 7.5 -e 0.001 --threads 10

#--min-seq-id (sequence identity cutoff)

#-c (coverage cutoff)

#-s (sensitivity, usually just stay at the max of 7.5)

#-e (evalue)