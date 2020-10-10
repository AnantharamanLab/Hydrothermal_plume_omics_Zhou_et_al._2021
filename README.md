# Hydrothermal_plume_omics_Zhou_et_al._2020

This GitHub repository contains bioinformatic analyzing methods within the manuscript XXX need to change_Ecology of inorganic sulfur auxiliary metabolism in widespread bacteriophages_

The bioRxiv preprint of this manuscript can be found [here](https://www.biorxiv.org/content/10.1101/2020.08.24.253096v1). 

10/6/2020  
Zhichao Zhou   
Karthik Anantharaman  
University of Wisconsin-Madison 

[Anantharaman Lab](https://anantharamanlab.com/)


## Table of Contents:
1. [Data files](#data)
2. [Phage analysis](#phage analysis)
3. [Phage to bacterial sulfur-related AMG ratio calculation](#map)
4. [Metatranscriptomic mapping](#ratio)
5. [Contact](#contact)


## Explanations


### Data Files <a name="data"></a>

* `Kieft_and_Zhou_et_al_2020.genomes.fasta`: all 191 vMAG genome sequences used in this study  
* `Kieft_and_Zhou_et_al_2020.proteins.faa`: all predicted protein sequences for the vMAGs used in this study  


### Phage analysis <a name="phage analysis"></a>
_Dependencies_: Perl v5.22.1, MMseqs2 Version: 7.4e23d, Python 3.5.2

(Higher version of each software should be OK, normally. Detailed information for backward compatibility could be found in the official website for each software )

_Explanation to each step in the pipeline:_  

**Step 1** 01.Get_statistics_from_VIBRANT_result.pl

Get the statistics of phage coding density from VIBRANT result. This will give a table with head line containing: "Plume ID", "Phage sequence number", "Metagenome sequence number", "Phage sequence density (n / 1,000 metagenome sequences)", "Phage nucleotide number", "Metagenome nucleotide number", "Phage nucleotide density (n / 1,000,000 metagenome nucleotides)". Resulted table indicates two phage coding densities (phage sequence density and phage nucleotide density) for each metagenome.

**Step 2** 02.Run_phage_taxonomy.pl and phage_taxonomy_tool_v6.py

Get phage taxonomy at subfamily level. Resulted folder contains two files: "~.VIVID.protein-taxonomy.tsv" and “~.VIVID.virus-taxonomy.tsv”, showing phage protein and phage genome taxonomy results respectively.

**Step 3** 03.Change_seq_head.pl

Change the head of sequences in a given fasta/faa file to 1) replace space to underscore, and cut head by its first "\t", and 2) add plume ID to sequence head. This helps to separate each metagenome, and make tracing each metagenome and each scaffold easier in the downstream analysis.

**Step 4** 04.Dereplicate_phage_proteins.sh

Dereplicate phage proteins first according to their sequence similarity.  The setting was "--min-seq-id 0.95 -c 0.8 -s 7.5 -e 0.001 ", which means sequence identity cutoff set as 95%, coverage cutoff set as 80%, sensitivity   coefficient set as 7.5 (usually just stay at the max of 7.5), and e-value set as 0.001.

**Step 5** 05.Run_mmseqs_to_all_phage_proteins.sh

Run MMseqs2 to cluster all phage proteins. The setting was "--min-seq-id 0.25 -c 0.5 -s 7.5 -e 0.001".

**Step 6** 06.Filter_small_clusters.sh

Filter small clusters which contain only proteins from 10 genomes. This will result a folder named "\*mmseqs_clusters", within which all cluster proteins were provided as "cluster_\*.faa".

**Step 7** 07.Make_genome_list.pl

Firstly, use "ls" to make a list for all cluster protein files in folder "\*mmseqs_clusters" and name the list as "cluster_list.txt". Second, use this script to make "genomes_list.txt" in folder "\*mmseqs_clusters".

**Step 8** 08.Count_cluster_hits_by_genome.sh and count_cluster_hits_by_genome.py

Count cluster hits by genome, and write result into a table with in folder "\*mmseqs_clusters".

**Step 9** 09.Make_metadata.pl

Make metadata table for all phage genomes. It contains columns "Genome", "Taxonomy", "Environment",	"Lysogenic property", "Clusters".

**Step 10** heatmap_plot.R

The R script needs open two input files and produce the heatmap to cluster phage genomes.  The "Count_cluster_hits_by_genome_result.tsv" generated by Step 8 and "Gene_and_taxonomy_metadata.tsv" generated by Step 9 will be used as inputs. The heatmap contains strips of metadata input, dendrogram of genome clusters, and dendrogram of protein clusters.

**Step 11** 10.Get_unknown_phage_proteins_from_each_metagenome.pl

Since the result of clustering all phage genomes are not good (~10k genomes), we try to only analyze the unknown phages (~5k genomes). This script helps to get unknown phage proteins from each metagenome and put them into a new folder.

**Step 12** 11.Dereplicate_phage_proteins.batch_for_each_metagenome.pl

Dereplicate the unknown phage proteins for each metagenome. Settings are the same with those in Step 4.

**Step 13** 12.Run_mmseqs_to_all_phage_proteins.batch_for_each_metagenome.pl

Run MMseqs2 to cluster all unknown phage proteins for each metagenome.  Settings are the same with those in Step 5.

**Step 14** 13.Summarize_protein_diversity.pl

Summarize the unknown phage protein diversity information. Produce a table for each metagenome containing the information of "Unknown phage protein num", "dRep95 unknown phage protein num", "Unknown phage protein cluster num".

The result suggests unknown phages share limited sequence similarity and can hardly be divided into certain groups by clustering phage proteins. Here, we only give the summary information to show the diverse phage protein pattern of unknown phages.



### Phage to bacterial sulfur-related AMG ratio calculation <a name="map"></a>

_Dependencies_: Perl v5+, Bowtie 2 v2.3.4.1,  jgi_summarize_bam_contig_depths (implemented in MetaWrap), Diamond v0.9.28.129, MAFFT v7.271, IQ-TREE v1.6.9

_Explanation to each step in the pipeline:_    
01.metagenome_mapping.pl    
Concatenate all the genes from metagenome as the mapping reference; Use Bowtie 2 to map filterd and QC-processed reads; Finally get sorted bam files as the result. The metadata file "Metagenome_map.txt" was used to allow processing mutiple metagenomes mapping by this script.    

02.calculate_metagenome_to_MAG_depth.pl    
Use "jgi_summarize_bam_contig_depths" to calculate gene coverage.    

03.grep_dsrA_list.sh    
Grep _dsrA_ genes from metagenome. Each IMG metagenome has its annotation by the DOE IMG database. We used the annotation to pre-select _dsrA_ genes from each metagenome. It needs further manual curation.  

04.calculate_viral_to_bacterial_sulfur-related_AMG_ratio.pl    
Read gene coverage result from Step #2. Calculate the viral to bacterial total sulfur-related AMG gene coverage ratio for each metagenome.    

05.grep_all_dsrA_gene.pl    
Grep all the dsrA gene encoding proteins from the metagenome.    

06.run_blastp.sh    
Run Diamond Blastp for all the DsrA sequences.    

07.parse_blast_result.sh    
Parse the Diamond Blastp result.    

08.find_top_non_virus_hit.pl    
Parse the Diamond Blastp result to screen for top 10 non-virus hits, which will be downloaded and used as reference sequences to build phylogenetic tree.    

09.make_dsrA_tree.pl    
Use viral and bacterial DsrA sequences, and reference DsrA sequences from Step #8, to build phylogenetic tree. MAFFT was used to align the sequences, and IQ-TREE was used to build the tree.    

10.Replace_tip_names.pl    
Replace tip names of resulted tree to the ones that are meaningful and formal.       

11.calculate_viral_to_bacterial_sulfur-related_AMG_ratio_according_to_taxonomy.pl
Firstly, divide the viral and bacterial sequences to each category according to their taxonomy. Then, calculate their AMG gene coverage ratios.    

12.calculate_viral_to_putative_host_bacterial_pair_sulfur-related_AMG_ratio.pl        
Firstly, get the viral to putative host bacterial pair information. Then, calculate viral and bacterial sequence gene coverage percentage values within each pair (the gene coverage values of viral and bacterial sequence were all normalized by gene numbers). Fianlly, calcualte the viral to putative host bacterial gene coverage ratios.


### Metatranscriptomic mapping <a name="ratio"></a>

_Dependencies_: Perl v5+, Bowtie 2 v2.3.4.1,  pileup (implemented in BBmap), Diamond v0.9.28.129, MAFFT v7.271, IQ-TREE v1.6.9    

01.Transcriptom_mapping.pl    
Concatenate all the genes from metagenome as the mapping reference; Use Bowtie 2 to map filterd, QC-processed and rRNA removed metatranscriptomic reads; Finally get sorted bam files as the result. The metadata file "Transcriptome_map.txt" was used to allow processing mutiple metatranscriptomes mapping by this script.        

02.pileup_to_calculate_gene_reads_abundance_for_MetaT.sh
Use "pileup.sh" to calculate gene abundance.   

03.parse_pileup_info.pl    
Parse pileup result.    

04.calculate_MetaT_RPKM.pl    
Calculate the gene expression result by normalizing gene abundance by read numbers (1 million reads) and gene length (1k).    

05.parse_MetaT_RPKM_result.pl    
Parse to get MetaT RPKM results for targeted genes.    

06.run_blastp.sh    
Run Diamond Blastp for all the DsrA sequences.    

07.parse_blast_result.sh    
Parse the Diamond Blastp result.    

08.find_top_non_virus_hit.pl     
Parse the Diamond Blastp result to screen for top 10 non-virus hits, which will be downloaded and used as reference sequences to build phylogenetic tree.     

09.make_dsrA_tree.pl    
Use viral and bacterial DsrA sequences, and reference DsrA sequences from Step #8, to build phylogenetic tree. MAFFT was used to align the sequences, and IQ-TREE was used to build the tree.     

10.Replace_tip_names.pl    
Replace tip names of resulted tree to the ones that are meaningful and formal.    

### Contact <a name="contact"></a>

Zhichao Zhou,  zczhou2017@gmail.com & zzhou388@wisc.edu
