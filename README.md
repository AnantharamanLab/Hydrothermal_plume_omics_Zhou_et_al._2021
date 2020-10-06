# Kieft_and_Zhou_et_al._2020

This GitHub repository contains bioinformatic analyzing methods within the manuscript XXX need to change_Ecology of inorganic sulfur auxiliary metabolism in widespread bacteriophages_

The bioRxiv preprint of this manuscript can be found [here](https://www.biorxiv.org/content/10.1101/2020.08.24.253096v1). 

10/6/2020  
Zhichao Zhou   
Karthik Anantharaman  
University of Wisconsin-Madison 

[Anantharaman Lab](https://https://anantharamanlab.com/)


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
_Dependencies_: Perl v5+, Bowtie 2 v2.3.4.1,  jgi_summarize_bam_contig_depths (implemented in MetaWrap), Diamond v0.9.28.129, MAFFT v7.271, IQ-TREE v1.6.9

_Explanation to each step in the pipeline:_  

01.Get_statistics_from_VIBRANT_result.pl

Get the statistics of phage coding density from VIBRANT result. This will give a table with head line containing: "Plume ID", "Phage sequence number", "Metagenome sequence number", "Phage sequence density (n / 1,000 metagenome sequences)", "Phage nucleotide number", "Metagenome nucleotide number", "Phage nucleotide density (n / 1,000,000 metagenome nucleotides)". Resulted table indicates two phage coding densities (phage sequence density and phage nucleotide density) for each metagenome.

02.Run_phage_taxonomy.pl and phage_taxonomy_tool_v6.py

Get phage taxonomy at subfamily level. Resulted folder contains two files: "~.VIVID.protein-taxonomy.tsv" and “~.VIVID.virus-taxonomy.tsv”, showing phage protein and phage genome taxonomy results respectively.



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

Kristopher Kieft, kieft@wisc.edu  
Zhichao Zhou, zzhou388@wisc.edu 
