#!/bin/bash
#
# pipeline for RNA-seq data analysis:
# calls for subscripts.
##
##
#
#
START=$(date +%s)
#       check for installed programs : FastQC, TOPHAT, Cutadapt, STAR
#
#
# http://stackoverflow.com/questions/14637979/how-to-permanently-set-path-on-linux
# echo "export PATH=$PATH:/home/talwar/bin" >> ~/.profile
# source ~/.profile 
# ln -s ~/installed_applications/STAR-STAR_2.4.1d/source/STAR ~/bin/
#
#
#
if ! which fastqc > /dev/null || ! which cutadapt > /dev/null || ! which tophat2 > /dev/null || ! which STAR > /dev/null || ! which htseq-count > /dev/null || ! which featureCounts > /dev/null
then
echo -e "    Either ONE or ALL of the programmes not available: 'fastqc', 'cutadapt', 'tophat', 'star', 'HTSeq-counts', 'FeatureCounts' please check if all the required programmes are installed and run the script again!
    quitting...."
exit
else
echo " "
echo " " 
echo "                                                                                `tput bold`All required programmes available!"
echo "   " 
#
#
#
echo "    `tput bold`The analysis will be performed in this sequence:
    0. creating directory structure (Automatic)
    1. Quality check (fastqc)
    2. Trimmming data using cutadapt
    3. Mapping data using TopHat
    4. Mapping data using STAR
    5. counting data using HT-Seq & FeatureCount
    6. differential expression analysis using DESeq"
#
#
#
#
echo " "
echo "Creating Directory structure :" 
mkdir -p ~/pool-talwar/RNA-Seq_pipeline
mkdir -p ~/pool-talwar/RNA-Seq_pipeline/Quality_control
mkdir -p ~/pool-talwar/RNA-Seq_pipeline/Quality_control/fastqc
mkdir -p ~/pool-talwar/RNA-Seq_pipeline/Quality_trimmed_files
mkdir -p ~/pool-talwar/RNA-Seq_pipeline/Quality_trimmed_files/Pair1
mkdir -p ~/pool-talwar/RNA-Seq_pipeline/Quality_trimmed_files/Pair2
mkdir -p ~/pool-talwar/RNA-Seq_pipeline/Raw_data_files
mkdir -p ~/pool-talwar/RNA-Seq_pipeline/Raw_data_files/pair1
mkdir -p ~/pool-talwar/RNA-Seq_pipeline/Raw_data_files/pair2
mkdir -p ~/pool-talwar/RNA-Seq_pipeline/scripts
mkdir -p ~/pool-talwar/RNA-Seq_pipeline/Mapping
mkdir -p ~/pool-talwar/RNA-Seq_pipeline/Mapping/HT-Seq_counts
mkdir -p ~/pool-talwar/RNA-Seq_pipeline/Mapping/featurecounts
mkdir -p ~/pool-talwar/RNA-Seq_pipeline/Differential_Expression
#
#
#
echo " "
echo " "
#
echo "   ***** Important : Now place the raw data files in folders : ~/pool-talwar/RNA-Seq_pipeline/Raw_data_files/pair1 & pair 2 
    for example : place the read pair 1 in pair1 & read pair 2 in pair2 folder resp. " 
#
#
echo " " 
#
echo "    Press 1 followed by [Enter] when your done copying the raw files to the folders, to start Fastqc: "
##
#
##################################################################################################
#########################                                               ##########################
#########################              QUALITY-CHECK (FASTQC)           ##########################
##################################################################################################
#
##
echo "  "
read Input
if [ $Input -eq "1" ]
then
echo "    starting quality check at :`date`"
echo " "
sh ~/pool-talwar/RNA-Seq_pipeline/scripts/Quality_check.sh #>> ~/pool-talwar/RNA-Seq_pipeline/Readme_output.txt  #### the location for the quality control file.
echo "    `tput bold`Quality check (fastqc) ended at :`date``tput sgr0`"
echo "    "
echo "    `tput bold`XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX-------------------------------XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX`tput sgr0`"
echo "    "
else
echo "    sorry, you selected $Input, this is not the right option. Please press 1 and [Enter]`tput sgr0`"
#
fi
echo " "
echo "     `tput bold` Create a file named 'quality_trimming_parameters.txt' with first column as quality + ' ' + cutoff + ' ' + value(for example 28) 
     and place it in the directory : ~/home/talwar/pool-talwar/RNA-seq_pipeline/"
echo "  " 
echo "    `tput bold`Press 2 followed by [Enter] to start QUALITY-TRIMMING (CUTADAPT): "
read Input
if [ $Input -eq "2" ]
then
#
#######################################################################################################
#########################                                                    ##########################
#########################              QUALITY-TRIMMING (CUTADAPT)           ##########################
#######################################################################################################
##
echo "    `tput bold`Starting quality trimming now!"
echo "    starting quality trimming at :`date``tput sgr0`"
echo " "
sh ~/pool-talwar/RNA-Seq_pipeline/scripts/trimming.sh  #>> ~/pool-talwar/RNA-Seq_pipeline/Readme_output.txt #### the location for the quality trimming file.
echo "    `tput bold`Quality trimming (cutadapt) ended at :`date`"
echo "    "
echo "    `tput bold`XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX-------------------------------XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX`tput sgr0`"
echo "    "
##
else
echo "    sorry, you selected $Input, this is not the right option. Please press 2 and [Enter]`tput sgr0`"
#
fi
echo " "
echo " "
echo "     `tput bold` Create a file named 'mapping_parameters.txt' with first column as maximum + ' ' + multihits + ' ' + value(for example 3) 
     and place it in the directory : ~/home/talwar/pool-talwar/RNA-seq_pipeline/"
echo "  " 
#
#######################################################################################################
#########################                                                    ##########################
#########################              MAPPING (TOPHAT)                      ##########################
#######################################################################################################
##
echo " "
echo "    `tput bold`Type 3 followed by [Enter] to start MAPPING (TOPHAT) or type 4 followed by [Enter] to start MAPPING (STAR): "
read Input
if [ $Input -eq "3" ]
then
echo "    `tput bold`"
echo "    starting mapping (TopHat) at :`date`"
echo " "
sh ~/pool-talwar/RNA-Seq_pipeline/scripts/mapping.sh #>> ~/pool-talwar/RNA-Seq_pipeline/Readme_output.txt  #### the location for the mapping file.
echo "    `tput bold`Mapping (TopHat) ended at :`date``tput sgr0`"
echo "    "
echo "    `tput bold`XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX-------------------------------XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX`tput sgr0`"
#
##
#######################################################################################################
#########################                                                    ##########################
#########################              MAPPING (STAR)                        ##########################
#######################################################################################################
##
#
elif [ $Input -eq "4" ]
then
echo "    `tput bold`"
echo "    starting mapping (STAR) at :`date`"
echo " "
sh ~/pool-talwar/RNA-Seq_pipeline/scripts/mapping_STAR.sh #>> ~/pool-talwar/RNA-Seq_pipeline/Readme_output.txt  #### the location for the mapping file.
echo "    `tput bold`Mapping (STAR) ended at :`date``tput sgr0`"
echo "    "
echo "    `tput bold`XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX-------------------------------XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX`tput sgr0`"
echo "    "
#
#
else
echo "    sorry, you selected $Input, this is not the right option. Please press 2 or 3 and [Enter]`tput sgr0`"
#
fi
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "    `tput bold`It took $DIFF seconds to finish the whole run! the output readme file with all outputs can be found with name : Readme_output.txt`tput sgr0`"
#
#
# finish below is for the if statement line #20
echo " "
fi
#
##
#
