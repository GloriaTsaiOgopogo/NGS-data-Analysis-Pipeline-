#!/bin/bash
#
#
#### Mapping script for paired end data 
#
[ -f ~/pool-talwar/RNA-Seq_pipeline/scripts/mapping_parameters.txt ] && echo "    `tput bold`parameters file found!`tput sgr0`" || echo "    `tput bold`parameter file Not found! Please create a file with name: mapping_parameters.txt`tput sgr0`"
#
#    For Mapping :
echo " "
echo "    `tput bold`Mapping with Tophat will now begin :`tput sgr0`"
for file1 in ~/pool-talwar/RNA-Seq_pipeline/Quality_trimmed_files/Pair1/*.fastq
do
while read f1 f2 f3
do
f3=$f3
PAIR1=$(echo $file1 | sed -E 's/\/home\/talwar\/pool-talwar\/RNA-Seq\_pipeline\/Quality\_trimmed\_files\/Pair1\/(.*)_1_trimmed.fastq/\1/')
PAIR2=$(echo $file1 | sed -E 's/\/home\/talwar\/pool-talwar\/RNA-Seq\_pipeline\/Quality\_trimmed\_files\/Pair1\/(.*)\_1_trimmed/\1\_2_trimmed/g')
#
echo " "
#echo "    $PAIR1"
#echo "    $PAIR2"
echo "    `tput bold`Maximum multihits value is: $f3`tput sgr0`"
echo "    `tput bold`Running Tophat for $PAIR1`tput sgr0`"
echo " "
echo "tophat2 -o ~/pool-talwar/RNA-Seq_pipeline/Mapping/$PAIR1\_tophat_out/ -p 13 -g $f3 -G ~/Masters_Thesis_Project/Genomes/Mus_musculus/GRCm38/Sequence/Bowtie2Index/genome.gtf ~/Masters_Thesis_Project/Genomes/Mus_musculus/GRCm38/Sequence/Bowtie2Index/genome $file1 ~/pool-talwar/RNA-Seq_pipeline/Quality_trimmed_files/Pair2/$PAIR2"
done < ~/pool-talwar/RNA-Seq_pipeline/scripts/mapping_parameters.txt
#
#                                             ###############################
#                                             ###     Running HT-seq      ###
#                                             ###############################
#
# insert parameter file for HTseq specifying the ID featuretype
#
echo " "
echo "                                             `tput bold`Running HTSeq-counts now `tput sgr0`"
echo " "
#mv ~/pool-talwar/RNA-Seq_pipeline/Mapping/$PAIR1\_tophat_out/accepted_hits.bam ~/pool-talwar/RNA-Seq_pipeline/Mapping/$PAIR1\_tophat_out/$PAIR1\_accepted_hits.bam
#
# Sort the bam files and convert them to sam files.
#
#samtools sort -n -@ 13 ~/pool-talwar/RNA-Seq_pipeline/Mapping/$PAIR1\_tophat_out/$PAIR1\_accepted_hits.bam ~/pool-talwar/RNA-Seq_pipeline/Mapping/$PAIR1\_tophat_out/$PAIR1\_accepted_hits
#
#
#samtools index ~/pool-talwar/RNA-Seq_pipeline/Mapping/$PAIR1\_tophat_out/$PAIR1\_accepted_hits.bam
#
# converting to sam for input for htseq-count
#
#samtools view ~/pool-talwar/RNA-Seq_pipeline/Mapping/$PAIR1\_tophat_out/$PAIR1\_accepted_hits.bam >  ~/pool-talwar/RNA-Seq_pipeline/Mapping/$PAIR1\_tophat_out/$PAIR1\_accepted_hits.sam
#
##
#
#htseq-count ~/pool-talwar/RNA-Seq_pipeline/Mapping/$PAIR1\_tophat_out/$PAIR1\_accepted_hits.sam ~/Masters_Thesis_Project/Genomes/Mus_musculus/GRCm38/Sequence/Bowtie2Index/genome.gtf > ~/pool-talwar/RNA-Seq_pipeline/Mapping/$PAIR1\_tophat_out/$PAIR1\_count.txt
#
### copy the files to another location as well:
#
#cp ~/pool-talwar/RNA-Seq_pipeline/Mapping/$PAIR1\_tophat_out/$PAIR1\_count.txt ~/pool-talwar/RNA-Seq_pipeline/Mapping/HT-Seq_counts/
#
#####                    Rename the accepted_hits.bam files
#
#
#ln -s ~/pool-talwar/RNA-Seq_pipeline/Mapping/$PAIR1\_tophat_out/$PAIR1\_accepted_hits.bam ~/pool-talwar/RNA-Seq_pipeline/Mapping/featurecounts/
#rm ~/pool-talwar/RNA-Seq_pipeline/Mapping/$PAIR1\_tophat_out/$PAIR1\_accepted_hits.sam
done
echo " " 
echo "     HT-Seq count is now complete! files are at: ~/pool-talwar/RNA-Seq_pipeline/Mapping/HT-Seq_counts/"
#
#                                             ###############################
#                                             ###  Running Feature-count  ###
#                                             ###############################
#
echo " "
echo "                                             `tput bold`Running Feature-count now `tput sgr0`"
echo " "
#feature_count_input=`ls ~/pool-talwar/RNA-Seq_pipeline/Mapping/featurecounts/`
echo "    taking these files for feature count as input : 
    $feature_count_input"
#echo "$feature_count_input" | tr '\n' ' '
#
#        Running featurecounts:
#
echo "featureCounts -a ~/Masters_Thesis_Project/Genomes/Mus_musculus/GRCm38/Sequence/Bowtie2Index/genome.gtf -o ~/pool-talwar/RNA-Seq_pipeline/Mapping/featurecounts/Complete_featurecounts.txt -T 13 -P -p -d 25 -D 101 -B -C ~/pool-talwar/RNA-Seq_pipeline/Mapping/featurecounts/*.bam"
#
#
# cp ~/pool-talwar/RNA-Seq_pipeline/Mapping/HT-Seq_counts/*.txt ~/pool-talwar/RNA-Seq_pipeline/Differential_Expression/
# cp ~/pool-talwar/RNA-Seq_pipeline/Mapping/featurecounts/Complete_featurecounts.txt ~/pool-talwar/RNA-Seq_pipeline/Differential_Expression/
echo " " 
echo "    `tput bold`Mapping is now complete!, the output files can be found at : ~/pool-talwar/RNA-Seq_pipeline/Mapping/"
echo " "
echo "     `tput bold`Feature count is also complete! files are at: ~/pool-talwar/RNA-Seq_pipeline/Mapping/featurecounts/`tput sgr0`"
echo " "
#
#
#
#################                XXXXXX-------------------------XXXXXXXXXX                  ##########################
