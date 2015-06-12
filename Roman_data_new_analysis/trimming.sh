### script for trimming data after quality check.
#
#
#                                           First to check for presence of adapter in each file. 
#                                         searching is only done for illumina adapter 'AGATCGGAAGAGC'
#
echo "   `tput bold`checking for illumina adapter('AGATCGGAAGAGC') in : PAIR1`tput sgr0`"
echo " "
for file1 in /home/talwar/pool-talwar/RNA-Seq_pipeline/Raw_data_files/pair1/*.fastq
do
PAIR1=$(echo $file1 | sed -E 's/\/home\/talwar\/pool-talwar\/RNA-Seq\_pipeline\/Raw\_data\_files\/pair1\/(.*)/\1/')
echo "    `tput bold`$PAIR1`tput sgr0`"
#grep --colour=always -o 'AGATCGGAAGAGC' $file1 | wc -l
done
#
###   check for illumina adapter in piar2 :
#
echo "    `tput bold`checking for illumina adapter('AGATCGGAAGAGC') in : PAIR2`tput sgr0`"
echo " "
for file2 in ~/pool-talwar/RNA-Seq_pipeline/Raw_data_files/pair2/*.fastq
do
PAIR2=$(echo $file2 | sed -E 's/\/home\/talwar\/pool-talwar\/RNA-Seq\_pipeline\/Raw\_data\_files\/pair2\/(.*)/\1/')
echo "    `tput bold`$PAIR2`tput sgr0`"
#grep --colour=always -o 'AGATCGGAAGAGC' $file2 | wc -l
done
#
#
#                                            ##############################################
#                                            ##  Now for quality trimming using CUTADAPT ## 
#                                            ##############################################
#
#
echo " "
echo "    `tput bold`quality trimming will now begin!`tput sgr0`"
echo " "
echo "    `tput bold`checking for parameter file now`tput sgr0`"



#######################################
### if you wish to use user input then use these commands : 
#echo "    `tput bold`Input the quality cutoff for MATE 1: `tput sgr0`"
#read parameter1
#if [ ! $parameter1 -eq "0" ]
#then
#param1=$parameter1
#else
#echo "     The file 'quality_trimming_parameters.txt' not found !"
###################################


##  here check for a file if it exists! which is named "quality_trimming_parameters.txt" 
##  if the files exists then read the parameters from each line !    specify the format of the file !

[ -f ~/pool-talwar/RNA-Seq_pipeline/scripts/quality_trimming_parameters.txt ] && echo "    `tput bold`parameters file found!`tput sgr0`" || echo "    parameter file Not found!"
#
#    For Mate 1 :
echo " "
echo " "
echo "    `tput bold`Running Cutadapt for MATE1:`tput sgr0`"
for file1 in /home/talwar/pool-talwar/RNA-Seq_pipeline/Raw_data_files/pair1/*.fastq
do
while read f1 f2 f3
do 
echo "    `tput bold`Quality cutoff value is: $f3`tput sgr0`"
echo " "
f3=$f3
PAIR1=$(echo $file1 | sed -E 's/\/home\/talwar\/pool-talwar\/RNA-Seq\_pipeline\/Raw\_data\_files\/pair1\/(.*)/\1/')
echo "    ~/bin/cutadapt --quality-cutoff $f3 --quiet -o ~/pool-talwar/RNA-Seq_pipeline/Quality_trimmed_files/$PAIR1\_trimmed.fastq $file1"
done < quality_trimming_parameters.txt
done
#
#
#
#
#    For MATE2:
echo "    `tput bold`Running Cutadapt for MATE2:`tput sgr0`"
for file2 in ~/pool-talwar/RNA-Seq_pipeline/Raw_data_files/pair2/*.fastq
do
PAIR2=$(echo $file2 | sed -E 's/\/home\/talwar\/pool-talwar\/RNA-Seq\_pipeline\/Raw\_data\_files\/pair2\/(.*)/\1/')
while read f1 f2 f3
do 
f3=$f3
echo "    ~/bin/cutadapt --quality-cutoff $f3 --quiet -o ~/pool-talwar/RNA-Seq_pipeline/Quality_trimmed_files/$PAIR2\_trimmed.fastq $file2"
done < quality_trimming_parameters.txt
done
#
#
echo " "
echo "    `tput bold`Quality trimming has ended !"
echo "    The output files can be found at : ~/pool-talwar/RNA-Seq_pipeline/Quality_trimmed_files/"
echo " "
tput sgr0
#
#
#
#