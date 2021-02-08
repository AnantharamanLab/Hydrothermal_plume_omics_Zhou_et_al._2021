# Hydrothermal_plume_omics_Zhou_et_al._2021

This GitHub repository contains bioinformatic analyzing methods within the manuscript:

XXX

The bioRxiv preprint of this manuscript can be found [here](https://www.biorxiv.org/content/10.1101/2020.08.24.253096v1). 



10/6/2020  
Zhichao Zhou   
Karthik Anantharaman  
Department of Bacteriology, University of Wisconsin-Madison 

[Anantharaman Lab link](https://anantharamanlab.com/)


## Table of Contents:

[Phage analysis](#phage_analysis)

[Metagenome and metatranscriptome mapping](#mapping)

[Comparison among three ecosystems](#comparison)

[Comparison within Mid-Cayman Rise ecosystem](#comparison_cayman)

[Key microbiome analysis input fasta files](#key_microbiome)

[Energy distribution calculation based on MetaGs](#energy_distribution)

[Contact](#contact)


## Explanations


### Phage analysis <a name="phage_analysis"></a>
_Dependencies_: Perl v5.22.1, MMseqs2 Version: 7.4e23d, Python 3.5.2

(Higher version of each software should be OK, normally. Detailed information for backward compatibility could be found in the official website for each software)

_Explanation to each step in the pipeline:_  

**Step 1** 01.Get_statistics_from_VIBRANT_result.pl

Get the statistics of phage coding density from VIBRANT result. This will give a table with head line containing: "Plume ID", "Phage sequence number", "Metagenome sequence number", "Phage sequence density (n / 1,000 metagenome sequences)", "Phage nucleotide number", "Metagenome nucleotide number", "Phage nucleotide density (n / 1,000,000 metagenome nucleotides)". Resulted table indicates two phage coding densities (phage sequence density and phage nucleotide density) for each metagenome.

**Step 2** 02.Run_phage_taxonomy.pl and phage_taxonomy_tool_v6.py

Get phage taxonomy at subfamily level. Resulted folder contains two files: "\*.VIVID.protein-taxonomy.tsv" and “*.VIVID.virus-taxonomy.tsv”, showing phage protein and phage genome taxonomy results respectively.

**Step 3** 03.Change_seq_head.pl

Change the head of sequences in a given fasta/faa file to 1) replace space to underscore, and cut head by its first "\t", and 2) add plume ID to sequence head. This helps to separate each metagenome, and make tracing each metagenome and each scaffold easier in the downstream analysis.

**Step 4** 04.Dereplicate_phage_proteins.sh

Dereplicate phage proteins first according to their sequence similarity.  The settings are "--min-seq-id 0.95 -c 0.8 -s 7.5 -e 0.001 ", which means sequence identity cutoff set as 95%, coverage cutoff set as 80%, sensitivity coefficient set as 7.5 (usually just stay at the max of 7.5), and e-value set as 0.001.

**Step 5** 05.Run_mmseqs_to_all_phage_proteins.sh

Run MMseqs2 to cluster all phage proteins. The settings are "--min-seq-id 0.25 -c 0.5 -s 7.5 -e 0.001".

**Step 6** 06.Filter_small_clusters_and_find_core_virome_proteins.sh and mmseqs-parse_core-unique-flexible.py

Filter small clusters which contain < 5 proteins and find core, flexible, and unique virome proteins.  This will result a resulted folder containing three subfolders. The minimum number of samples/locations for core virome proteins is 7. 

**Step 7** 07.Make_genome_list.pl

Firstly, use "ls" to make a list for all cluster protein files in core virome folder and name the list as "cluster_list.txt". Secondly, use this script to make "genomes_list.txt" in the folder. Since I combined Axial.Plume and Axial.Seawater into one Axial, so I modified the "genomes_list.txt" as "genomes_list.mdf.txt".

**Step 8** 08.Count_cluster_hits_by_genome.sh and count_cluster_hits_by_genome.py

Count core virome cluster hits by genome, and write result into a table within the folder

**Step 9** 09.Make_metadata.pl

Make metadata table for all phage genomes. It contains columns "Genome", "Taxonomy", "Environment",	"Lysogenic property", "Clusters".

**Step 10** heatmap_plot.R

The R script needs open two input files and produce the heatmap to cluster phage genomes.  The "Count_cluster_hits_by_genome_result.tsv" generated by Step 8 and "Gene_and_taxonomy_metadata.tsv" generated by Step 9 will be used as inputs. The heatmap contains strips of metadata input, dendrogram of genome clusters, and dendrogram of protein clusters.

**Step 11** 10.Get_unknown_phage_proteins_from_each_metagenome.pl

Since the result of clustering all phage genomes (~10k genomes) are not good , we try to only analyze the unknown phages (~5k genomes). This script helps to get unknown phage proteins from each metagenome and put them into a new folder.

**Step 12** 11.Dereplicate_phage_proteins.batch_for_each_metagenome.pl

Dereplicate the unknown phage proteins for each metagenome. Settings are the same with those in Step 4.

**Step 13** 12.Run_mmseqs_to_all_phage_proteins.batch_for_each_metagenome.pl

Run MMseqs2 to cluster all unknown phage proteins for each metagenome.  Settings are the same with those in Step 5.

**Step 14** 13.Summarize_protein_diversity.pl

Summarize the unknown phage protein diversity information. Produce a table for each metagenome containing the information of "Unknown phage protein num", "dRep95 unknown phage protein num", "Unknown phage protein cluster num".

The result suggests unknown phages share limited sequence similarity and can hardly be divided into certain groups by clustering phage proteins. Here, we only give the summary information to show the diverse phage protein pattern of unknown phages.



### Metagenome and metatranscriptome mapping <a name="mapping"></a>

Mapping metagenomic and metatranscriptomic reads separately for each environment

_Dependencies_: Perl v5.22.1, Bowtie 2 v2.3.4.1, jgi_summarize_bam_contig_depths (within metaWrap v1.0.2), pileup.sh (within BBmap)

(Higher version of each software should be OK, normally. Detailed information for backward compatibility could be found in the official website for each software)

_Explanation to each step in the pipeline:_   

**Step 1**  01.Transcriptom_mapping.Sep_Mapping.in_33.pl

Concatenate all the genes from genomes within each environment to make the mapping reference file. Use the QC-passed and rRNA-filtered metatranscriptomic reads to map against the reference file by Bowtie 2. The template input files, including "Transcriptom_map.mdf.txt" and "MAG_info.txt", are provided in the folder. In "Transcriptom_map.mdf.txt", metatranscriptomic reads are the reads that are labeled with "cDNA" in the second column.

**Step 2**  02.metagenome_mapping.Sep_Mapping.in_33.pl

Concatenate all the genomes within each environment to make the mapping reference file. Use the QC-passed metagenomic reads to map against the reference file by Bowtie 2. The template input files, including "Transcriptom_map.mdf.txt" and "MAG_info.txt", are provided in the folder.  In "Transcriptom_map.mdf.txt", metagenomic reads are the reads that are labeled with "DNA" in the second column.

**Step 3**  03.calculate_metagenome_to_MAG_depth.Sep_Mapping.pl

Calculate depth files based on "*.sorted.bam" files resulted from Step 2 (for metagenomes). "jgi_summarize_bam_contig_depths"  within metaWrap is used to do the calculation.

**Step 4**  04.pileup_to_calculate_gene_reads_abundance_for_MetaT.Sep_Mapping.sh

Calculate depth files based on "*.sorted.bam" files resulted from Step 1 (for metatranscriptomes). "pileup.sh" within BBmap is used to do the calculation.

**Step 5**  05.calculate_MAG_average_coverage.Sep_Mapping.pl

Calculate each MAG average coverage based on the result from Step 3. Resulted files are named as "*.MAG_average_coverage.txt". We normalize each metagenomic datasets to the size of 100M reads.

**Step 6**  06.parse_pileup_info.Sep_Mapping.pl

Parse the pileup outputs from Step 4.

**Step 7**  07.calculate_MetaT_TPM.Sep_Mapping.pl

Calculate each gene average MetaT expression level based on the result from Step 6. Resulted files are named as "*.MetaT.TPM.txt", and all the resulted values are in [TPM](https://rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/).



### Comparison among three ecosystems <a name="comparison"></a>

_Dependencies_: Perl v5.22.1, R version 4.0.2, R Studio Version 1.3.1073

(Higher version of each software should be OK, normally. Detailed information for backward compatibility could be found in the official website for each software)

_Explanation to each step in the pipeline:_    
**Step 1**  01.get_the_statistics_of_group_MetaG_abundance.pl

This script produces a table named "MetaG.10times.3_ecosystem.Group2MetaG_id_MetaG.xls" which contains the total gene coverage (10 times of the original gene coverage) of each microbial group. The gene coverage was used by multiplying 10, because the original numbers are too small. After this, we make a new table, named "MetaG.10times.3_ecosystem.Group2MetaG_id_MetaG.integer.txt" which transfers all the numbers from the previous table into integers. This will be used as the input for R script "DESeq2.MAG.coverage.3ecosystem.R ". The template input files, including "MAG_average_coverage.10times.txt" and "MAG_info.txt", are provided in the folder.

**Step 2**  02.get_the_statistics_of_group_MetaG_abundance_for_each_environment.pl

This script produces two tables named "MAG_average_coverage.Group2Row_Mean_MetaG.xls" and "MAG_average_coverage.Group2Row_Mean_MetaG_each_bin.xls". This script is the same as the Step 1 script, while it is adjusted for each environment.  The template input files, including "MAG_average_coverage.row_mean.simple.txt" and "MAG_info.txt", are provided in the folder for each environment. "MAG_average_coverage.row_mean.simple.txt" was made by taking the mean value of each row (including the background and plume mean values) for each environment based on "MAG_average_coverage.10times.txt". 

**Step 3**  03.parse_fun_normalized_abundance.Sep_Mapping.MetaT.v2.pl

Calculate normalized coverage for each functional trait. The input files for MAG metagenome coverage are from Step 5 of Metagenome and metatranscriptome mapping (values are transformed into percentages for each environment). Resulted files are named as "*.Fun2MetaG_abundance.txt". All the results are then combined into "Functional_analysis_summary.txt" (numbers are transferred into integers), which is provided in the folder.

**Step 4**  04.get_taxa_for_enriched_functions.pl

Get the microbial community contribution information to the enriched functions in each environment. The "enriched functions" (refer to "Enriched_functions.txt") are functions that are significantly enriched (having higher abundance) in each environment calculated by Step 3 and Rscript 2. The template input files, including "MAG_info.txt", "Fun_result.txt", "MetaG.10times_coverage.txt", and "Enriched_functions.txt", are provided in the folder. Resulted files are named as "Enriched_fun_micro_grp.*.txt".

**Step 5**  05.get_taxa_for_major_functions.pl

Get the microbial community contribution information to the major functions in each environment. The major functions include functions that are in the categories of carbon fixation, denitrification, sulfur cycling, hydrogen oxidation, and methane oxidation. The template input files, including "MAG_info.txt", "Fun_result.txt", "MetaG.10times_coverage.txt", and "Major_functions.txt", are provided in the folder. Resulted files are named as "Major_fun_micro_grp.*.txt".



*Rscripts*：

**Rscript 1** DESeq2.MAG.coverage.3ecosystem.R 
Perform DESeq2 analyses to compare the microbial community difference among different background and plume environment pairs, including B.Cym vs B.Lau, P.Cym vs P.Lau, P.Cym vs P.GyBn, and P.GyBn vs P.Lau. The Log2 Fold Change and *p*-value are calculated for each microbial group and included in the result.

**Rscript 2** DESeq2.MAG.function_analysis.3ecosystem.R 

Perform DESeq2 analyses to compare the microbial function difference among different background and plume environment pairs, including B.Cym vs B.Lau, P.Cym vs P.Lau, P.Cym vs P.GyBn, and P.GyBn vs P.Lau. The Log2 Fold Change and adjusted *p*-value are calculated for each microbial group and included in the result.

**Rscript 3** make_percentage_table.R

Make percentage table based on input table



### Comparison within Mid-Cayman Rise ecosystem <a name="comparison_cayman"></a>

_Dependencies_: Perl v5.22.1, Bowtie 2 v2.3.4.1, jgi_summarize_bam_contig_depths (within metaWrap v1.0.2), pileup.sh (within BBmap), R version 4.0.2, R Studio Version 1.3.1073

(Higher version of each software should be OK, normally. Detailed information for backward compatibility could be found in the official website for each software)

_Explanation to each step in the pipeline:_   

**Step 1**  01.Transcriptom_mapping_in_33.pl

Concatenate all the genes from genomes within Mid-Cayman Rise ecosystem (including Von Damm and Piccard) to make the mapping reference file (MAGs from Mid-Cayman Rise ecosystem were given in "Cayman_MAG.txt"). Use the QC-passed and rRNA-filtered metatranscriptomic reads to map against the reference file by Bowtie 2. The template input files, including "Transcriptom_map.txt" and "Cayman_MAG.txt", are provided in the folder. In "Transcriptom_map.txt", metatranscriptomic reads are the reads that are labeled with "cDNA" in the second column.

**Step 2**  02.metagenome_mapping_in_33.pl

Concatenate all the genomes within Mid-Cayman Rise ecosystem (including Von Damm and Piccard) to make the mapping reference file (The MAGs from Mid-Cayman Rise ecosystem were given in "Cayman_MAG.txt"). Use the QC-passed metagenomic reads to map against the reference file by Bowtie 2. The template input files, including "Transcriptom_map.txt" and "Cayman_MAG.txt", are provided in the folder.  In "Transcriptom_map.txt", metagenomic reads are the reads that are labeled with "DNA" in the second column.

**Step 3**  03.calculate_metagenome_to_MAG_depth.sh

Calculate depth files based on "*.sorted.bam" files resulted from Step 2 (for metagenomes). "jgi_summarize_bam_contig_depths" within metaWrap is used to do the calculation.

**Step 4**  04.pileup_to_calculate_gene_reads_abundance_for_MetaT.sh

Calculate depth files based on "*.sorted.bam" files resulted from Step 1 (for metatranscriptomes). "pileup.sh" within BBmap is used to do the calculation.

**Step 5**  05.calculate_MAG_average_coverage.pl

Calculate each MAG average coverage based on the result from Step 3. Resulted files are named as "*.MAG_average_coverage.txt". We normalize each metagenomic datasets to the size of 100M reads.

**Step 6**  06.parse_pileup_info.pl

Parse the pileup outputs from Step 4.

**Step 7**  07.calculate_MetaT_TPM.pl

Calculate each gene average MetaT expression level based on the result from Step 6. Resulted files are named as "*.MetaT.TPM.txt", and all the resulted values are in TPM.

**Step 8**  08.parse_fun_normalized_abundance_v2.pl

Calculate normalized coverage for each function trait based on metagenomes. The input files for MAG metagenome coverage are from Step 5 of metagenome mapping ("*.MAG_average_coverage.txt"). Resulted files are named as "Fun2MetaG_abundance.txt". 

**Step 9**  09.parse_fun_MetaT_TPM_v2.pl

Calculate normalized coverage for each function trait based on metatranscriptomes. The input files for MAG metatranscriptome coverage (or referred to as MAG expression level) files are from Step 7 of metatranscriptome mapping ("MetaT.TPM.txt"). Resulted files are named as "Fun2MetaT_abundance.txt". 

#### **Sub-analysis**:

1. ##### DESeq2 for *p* value

   **Step 1** 01.find_bin_info.pl

   Make bin list according to the input MAG information.

   **Step 2** 02.parse_MetaT_MAG_average_TPM.pl

   Parse the result of each gene average MetaT expression level ("*.MetaT.TPM.txt") into each MAG MetaT expression level ("Bin2MetaT_abundance.txt"). 

   **Step 3** 03.make.heatmap.table.01.CymD.MAG.MetaG.coverage.pl

   ​            03.make.heatmap.table.02.CymS.MAG.MetaG.coverage.pl

   ​            03.make.heatmap.table.03.CymD.MAG.MetaT.coverage.pl

   ​            03.make.heatmap.table.04.CymS.MAG.MetaT.coverage.pl

   ​            03.make.heatmap.table.08.CymS.Fun.MetaG.coverage.pl

   ​            03.make.heatmap.table.12.CymS.Fun.MetaT.coverage.pl

   (Some scripts are not given here, e.g., *CymD.Fun.MetaG* and *CymD.Fun.MetaT*, due to that there are no significant results from DESeq2 analyses)

   Make heatmap tables by using DESeq2-generated results (Rscript 1 results) which are used as inputs for "Heatmap.original.R".

   

2. ##### Taxa for enriched functions

   **Step 1** 01.get_taxa_for_enriched_functions.MetaG.pl

   Get the microbial community contribution information (based on metagenomes) to the enriched functions in each Mid-Cayman Rise sample. The "enriched functions" (refer to "Enriched_functions.txt") are functions that are significantly enriched (having higher abundance) in each Mid-Cayman Rise sample calculated by DESeq2 for *p* value Step 3. The template input files, including "MAG_info.Cayman.txt", "Fun_result.txt", "MAG_average_coverage.mdf.txt" (use the mean values of all background and plume samples, respectively), and "Enriched_functions.txt", are provided in the folder. Resulted files are named as "Enriched_fun_micro_grp.*.metaG.txt".

   **Step 2** 02.get_taxa_for_enriched_functions.MetaT.pl

   Get the microbial community contribution information (based on metatranscriptomes) to the enriched functions in each Mid-Cayman Rise sample. The "enriched functions" (refer to "Enriched_functions.metaT.txt") are functions that are significantly enriched (having higher abundance) in each Mid-Cayman Rise sample calculated by DESeq2 for *p* value Step 3. The template input files, including "MAG_info.Cayman.txt", "Fun_result.txt", "Bin2MetaT_abundance.mdf.txt" (use the mean values of all background and plume samples, respectively), and "Enriched_functions.metaT.txt", are provided in the folder. Resulted files are named as "Enriched_fun_micro_grp.*.metaT.txt".



*Rscripts*：

**Rscript 1** DESeq2.MAG.coverage.CymD.CymS.R

Perform DESeq2 analyses to compare the microbial community and function difference within Mid-Cayman Rise background and plume pairs, including: *CymD.MAG.MetaG/CymS.MAG.MetaG* (background and plume comparison for MAG composition based on metagenomes for CymD and CymS),  *CymD.MAG.MetaT/CymS.MAG.MetaT* (background and plume comparison for MAG composition based on metatranscriptomes for CymD and CymS),  *CymD.Fun.MetaG/CymS.Fun.MetaG* (background and plume comparison for functional trait based on metagenomes for CymD and CymS), and *CymD.Fun.MetaT/CymS.Fun.MetaT* (background and plume comparison for functional trait based on metatranscriptomes for CymD and CymS). Here, CymD stands for Piccard (Cayman Deep) and CymS stands for Von Damm (Cayman Shallow). The Log2 Fold Change and adjusted *p*-value are calculated for each microbial group and included in the result.

We provided the input files for this R script, which are "CymD.10times.transpose.MAG.MetaG.coverage.txt", "CymS.10times.transpose.MAG.MetaG.coverage.txt", "CymD.100times.MAG.MetaT.coverage.txt", "CymS.100times.MAG.MetaT.coverage.txt", "CymD.Fun2MetaG.txt", “CymS.Fun2MetaG.txt”, "CymD.Fun2MetaT.txt", and "CymS.Fun2MetaT.txt". These are modified subsets from the original files of "MAG_average_coverage.txt" (from Step 5), "Bin2MetaT_abundance.txt" (from Step7 and DESeq2 for *p* value Step 2), "Fun2MetaG_abundance.txt" (from Step 8), and "Fun2MetaT_abundance.txt" (from Step 9). Since that the numbers in input files should be integers, we multiplied the original MAG abundance and expression abundance by 10 times and 100 times as one can find in the name of input files, e.g., "*10times.transpose*" and "*100times*".  

**Rscript 2** Heatmap.original.R

Make heatmaps based on input files.



### Key microbiome analysis input fasta files <a name="key_microbiome"></a>

_Dependencies_:  QIIME v.1.9.1 

(QIIME was used here but not QIIME 2 due to that we only processed fasta files but not fastq files)

We put all the fasta files into the folder of "Key_microbiome_analysis_input_fasta_files".  One could recover and extract the tar.gz file like this:

cat fasta_files.tar.gz* > fasta_files.tar.gz ; tar xzf fasta_files.tar.gz

Routine [QIIME](http://qiime.org/) 16S rRNA gene analyzing method was used to interpret the microbial community distribution pattern. Please visit the official website for details.



### Energy distribution calculation based on MetaGs <a name="energy_distribution"></a>

_Dependencies_: Perl v5.22.1, Excel (Microsoft 365 Apps)

(Higher version of each software should be OK, normally. Detailed information for backward compatibility could be found in the official website for each software)

_Explanation to each step in the pipeline:_   

**Step 1** 01.parse_electron_transfer_MetaG_fun_norm_abund_v2.pl

Get the abundance of each electron transferring reaction based on metagenomes. The abundance of each electron transferring reaction was calculated by adding up coverage values of all genes that are responsible for the reaction. The abundance values were also normalized by gene numbers (aka. the enzyme numbers for each reaction). 

The input files include "*.MAG_average_coverage.txt" (store each MAG average coverage for different metagenomes; values were normalized by 100M reads) and "FunMap.list" (store the functional trait ID to each electron transferring reaction map). Template files are given for these input files.

**Step 2** Calculate energy contribution for each electron donor

Use Excel formula to calculate energy contribution for each electron donor. It was calculated by multiplying reaction abundance by the energy yield for each reaction. The reaction and energy yield information refer to Supplementary Table S2. The Excel template for calculating energy contribution is also placed in the folder ("Energy_contribution_calculation_template.xlsx"). The final curated excel file has been provided as one of the supplementary materials.

We used a similar calculation method for energy contribution based on metatranscriptomes. Detailed are not provided here.

​    

### Contact <a name="contact"></a>

Zhichao Zhou,  zczhou2017@gmail.com & zzhou388@wisc.edu

Karthik Anantharaman,   karthik@bact.wisc.edu 
