#!/bin/bash
## script for quality check for each file in the directory
## note that files should be in defined format if they are 'paired end' 
## the files should be unzipped and be in 'fastq' format
#
#
for file1 in ~/pool-talwar/RNA-Seq_pipeline/Raw_data_files/pair1/*.fastq
do
PAIR1=$file1
echo " "
echo "    `tput bold`Running fastqc for: $PAIR1`tput sgr0`"
echo "    fastqc -o ~/pool-talwar/RNA-Seq_pipeline/Quality_control/fastqc/ -t 13 -q $PAIR1"
done
#
#
###                              now for pair2:
#
#
for file2 in ~/pool-talwar/RNA-Seq_pipeline/Raw_data_files/pair2/*.fastq
do
PAIR2=$file2
echo " "
echo "    `tput bold`Running fastqc for: $PAIR2`tput sgr0`"
echo "     fastqc -o ~/pool-talwar/RNA-Seq_pipeline/Quality_control/fastqc/ -t 13 -q $PAIR2"
done
#
echo " "
echo "`tput bold`the FastQC outputs are at : ~/pool-talwar/RNA-Seq_pipeline/Quality_control/fastqc/ `tput sgr0`"
#
#
#
