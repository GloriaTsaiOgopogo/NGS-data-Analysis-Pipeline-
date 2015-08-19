#!/bin/bash
#
#### Mapping script for paired end data with STAR aligner 
#
[ -f ~/pool-talwar/RNA-Seq_pipeline/scripts/mapping_parameters.txt ] && echo "    `tput bold`parameters file found!`tput sgr0`" || echo "    `tput bold`parameter file Not found! Please create a file with name: mapping_parameters.txt`tput sgr0`"
#
#    For Mapping :
echo " "
echo "    `tput bold`Mapping with 'STAR' will now begin :`tput sgr0`"
for file1 in ~/pool-talwar/RNA-Seq_pipeline/Quality_trimmed_files/Pair1/*.fastq
do
while read f1 f2 f3
do
f3=$f3
PAIR1=$(echo $file1 | sed -E 's/\/home\/talwar\/pool-talwar\/RNA-Seq\_pipeline\/Quality\_trimmed\_files\/Pair1\/(.*)_1_trimmed.fastq/\1/')
PAIR2=$(echo $file1 | sed -E 's/\/home\/talwar\/pool-talwar\/RNA-Seq\_pipeline\/Quality\_trimmed\_files\/Pair1\/(.*)\_1_trimmed/\1\_2_trimmed/g')
#
echo " "
echo "    `tput bold`Maximum multihits value is: $f3`tput sgr0`"
echo "    `tput bold`Running STAR for $PAIR1, $PAIR2`tput sgr0`"
echo " "
echo "    $PAIR1"
#echo "    $PAIR2"
#mkdir ~/pool-talwar/RNA-Seq_pipeline/Mapping/$PAIR1\_STAR_out/
echo "    STAR --genomeDir ~/Masters_Thesis_Project/Genomes/Mus_musculus/STAR_genome/ --readFilesIn $file1 /home/talwar/pool-talwar/RNA-Seq_pipeline/Quality_trimmed_files/$PAIR2 --runThreadN 13 --outFileNamePrefix ~/pool-talwar/RNA-Seq_pipeline/Mapping/$PAIR1\_STAR_out/ --outReadsUnmapped Fastx"
done < ~/pool-talwar/RNA-Seq_pipeline/scripts/mapping_parameters.txt
#
#                                             ###############################
#                                             ###     Running HT-seq      ###
#                                             ###############################
#
echo " "
echo "                                                                                    `tput bold`Running HTSeq-counts now `tput sgr0`"
echo " "
echo "    htseq-count -q ~/pool-talwar/RNA-Seq_pipeline/Mapping/$PAIR1\_STAR_out/Aligned.out.sam ~/Masters_Thesis_Project/Genomes/Mus_musculus/GRCm38/Sequence/Bowtie2Index/genome.gtf > ~/pool-talwar/RNA-Seq_pipeline/Mapping/$PAIR1\_STAR_out/$PAIR1\_count.txt"
#
#
echo "    cp ~/pool-talwar/RNA-Seq_pipeline/Mapping/$PAIR1\_STAR_out/$PAIR1\_count.txt ~/pool-talwar/RNA-Seq_pipeline/Mapping/HT-Seq_counts/"
#
#####                    Rename the accepted_hits.bam files
echo "    mv ~/pool-talwar/RNA-Seq_pipeline/Mapping/$PAIR1\_STAR_out/Aligned.out.sam ~/pool-talwar/RNA-Seq_pipeline/Mapping/$PAIR1\_STAR_out/$PAIR1\_Aligned.out.sam"
echo "    ln -s ~/pool-talwar/RNA-Seq_pipeline/$PAIR1\_STAR_out/$PAIR1\_Aligned.out.sam ~/pool-talwar/RNA-Seq_pipeline/Mapping/featurecounts/"
done
echo "    `tput bold`Mapping is now complete!, the output files can be found at : ~/pool-talwar/RNA-Seq_pipeline/Mapping/"
echo "     HT-Seq count is also complete! files are at: ~/pool-talwar/RNA-Seq_pipeline/Mapping/HT-Seq_counts/`tput sgr0`"
echo " "
#
##                                             ###############################
#                                             ###  Running Feature-count  ###
#                                             ###############################
#
echo " "
echo "                                                                                     `tput bold`Running Feature-count now `tput sgr0`"
echo " "
#feature_count_input=`ls ~/pool-talwar/RNA-Seq_pipeline/Mapping/featurecounts/`
echo " taking these fiels for feature count as input : $feature_count_input"
#
#        Running featurecounts:
#
echo "    featureCounts -a ~/Masters_Thesis_Project/Genomes/Mus_musculus/mm10/Sequence/Bowtie2Index/genome.gtf -o ~/pool-talwar/RNA-Seq_pipeline/Mapping/featurecounts/Complete_featurecounts_STAR.txt -T 13 -P -C ~/pool-talwar/RNA-Seq_pipeline/Mapping/featurecounts/*.sam"
#
#
echo "     `tput bold`Feature count is also complete! files are at: ~/pool-talwar/RNA-Seq_pipeline/Mapping/featurecounts/`tput sgr0`"
echo " "
#
#
#################                XXXXXX-------------------------XXXXXXXXXX                  ##########################
