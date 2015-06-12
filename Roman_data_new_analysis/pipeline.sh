# pipeline for RNA-seq data analysis:
# calls for subscripts.

echo "`tput bold`    the analysis will be performed in this sequence:
   #creating directory structure= 1
    Quality check (fastqc) = 1
    Trimmming data using cutadapt = 2
    Mapping data using TopHat = 3
    Mapping data using STAR = 4 
    counting data = 5"
## getting user input for the step of analysis
##
##    Directory structure 
#mkdir -p ~/pool-talwar/RNA-Seq_pipeline
#mkdir -p ~/pool-talwar/RNA-Seq_pipeline/Quality_control
#mkdir -p ~/pool-talwar/RNA-Seq_pipeline/Quality_control/fastqc
#mkdir -p ~/pool-talwar/RNA-Seq_pipeline/Quality_control/trimming
#mkdir -p ~/pool-talwar/RNA-Seq_pipeline/Raw_data_files
#mkdir -p ~/pool-talwar/RNA-Seq_pipeline/Raw_data_files/pair1
#mkdir -p ~/pool-talwar/RNA-Seq_pipeline/Raw_data_files/pair2
#mkdir -p ~/pool-talwar/RNA-Seq_pipeline/scripts
#
echo "    Press 1 followed by [Enter] if you wihs to start the analysis for Quality check :"
##
#
##################################################################################################
#########################                                               ##########################
#########################              QUALITY-CHECK (FASTQC)           ##########################
##################################################################################################
#
##
read Input
if [ $Input -eq "1" ]
then
echo "    You selected Quality check 
"
echo "    starting quality check at :`date`"
sh ~/pool-talwar/RNA-Seq_pipeline/scripts/Quality_check.sh   #### the location for the quality control file.
echo "    Quality check (fastqc) ended at :`date`"
#tput sgr0
else
echo "    sorry, you selected $Input, this is not the right option. Please press 1 and [Enter]`tput sgr0`"
fi
##
#
#######################################################################################################
#########################                                                    ##########################
#########################              QUALITY-TRIMMING (CUTADAPT)           ##########################
#######################################################################################################
##
echo "    type in 2 followed by [enter] if you wish to perform quality trimming`tput sgr0`"
read Input
if [ $Input -eq "2" ]
then
echo "    starting quality trimming at :`date`"
sh ~/pool-talwar/RNA-Seq_pipeline/scripts/trimming.sh   #### the location for the quality trimming file.
echo "    Quality trimming (cutadapt) ended at :`date`"
tput sgr0
else
echo "    sorry, you selected $Input, this is not the right option. Please press 2 and [Enter]`tput sgr0`"
fi
##
#
#######################################################################################################
#########################                                                    ##########################
#########################              MAPPING (TOPHAT)                      ##########################
#######################################################################################################
##























