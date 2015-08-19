# RNA-Seq-Analysis-Pipeline-
scripts for automating the RNA-seq analysis of large datasets.

This script was generated using raw RNA-seq data from mouse samples after and before calories restriction.
the reads were paired end. Each mate had 50 M reads.
the steps that the script include for the analysis are :

1. FastQC (fasqc)
2. Adapter trimming and Quality trimming ( Cutadapt)
3. Mapping (Tophat2)
4. Read counts (HTseq)
5. differential expression analysis (DESeq2)


see the file ROman_data_new for the newest scripts for analysis of paired end data!

