#!/bin/sh

# multi_bwape_analysis_v1.sh
# 
#

# This script is designed to allow multiple samples/lanes of paired-end illumina data to be passed into the "BWApe_hg18_v1" pipeline

###############################################################################################################################
## To facilitate its use you must put uniquely named Illumina 1.3+ files in ngs/bwape/inputsequences/hold/lane(X)/read(1-2)   #
## It is essential that these file names are uniquely name or overwriting will occur										  #
## These files MUST have the ".txt" extension characteristic of the Illumina V1.3+ output "s_x_sequences.txt"				  #
###############################################################################################################################

#In this step we check that you are lauching the script from the correct location in case you are using it from a path directory

echo "`tput bold` ***Checking Current Directory is Correct***"

#List of directoryies to check
temp1=/export/Assa/ngs

#Checking if launch location is correct

if [ "`pwd`" != "$temp1" ] 
	then 
	echo " The script must be launched from the NGS directory "
	echo " The script was automatically killed due to a launch error - See Above Error Message" 
	exit 2                              
fi
echo ***Current	Directory is Correct***

#Check if files exist in the dilp hold folder

#List of directoryies to check
temp2=/export/Assa/ngs/bwape/inputsequences/hold/dilp/read1
temp3=/export/Assa/ngs/bwape/inputsequences/hold/dilp/read2

echo ***Checking dilp Hold Folder***
if [ `ls $temp2 | wc -l` != 1 ]       
	then 
	echo " The dilp Read1 hold folder does not contain the expect single file "
	echo " ERROR - The script was automatically killed due to a launch error - See Above Error Message"
	exit 2                     
fi
if [ `ls $temp3 | wc -l` != 1 ]       
	then 
	echo " The dilp Read2 hold folder does not contain the expect single file "
	echo " ERROR - The script was automatically killed due to a launch error - See Above Error Message"
	exit 2                     
fi
echo ***Found Expected Files***
#Current directory=ngs
echo ***Starting The Analysis of dilp***
date '+%m/%d/%y %H:%M:%S'
echo ***Starting The Analysis of dilp*** >> analysisnotes/Analysis.log
date '+%m/%d/%y %H:%M:%S' >> analysisnotes/Analysis.log

#In the next step we move the "dilp" data from ngs/bwape/inputsequences/hold/dilp/(read1-2) to ngs/bwape/inputsequences/illumina/(read1-2)

cd bwape/inputsequences/hold/dilp/read1
#Current Directory=ngs/bwape/inputsequences/hold/dilp/read1
echo Moving dilp Read1 File to Read1 Analysis Directory
date '+%m/%d/%y %H:%M:%S'
echo Moving dilp Read1 File to Read1 Analysis Directory >> ../../../../../analysisnotes/Analysis.log
date '+%m/%d/%y %H:%M:%S' >> ../../../../../analysisnotes/Analysis.log
echo Moving the following files:
for ligne in `ls *.txt`
do                                                                     
echo $ligne
done
for ligne in `ls *.txt`
do
echo Moving the following files: >> ../../../../../analysisnotes/Analysis.log
done
for ligne in `ls *.txt`
do
mv $ligne ../../../illumina/read1/
done
echo Moving dilp Read1 File to Read1 Analysis Directory Complete
date '+%m/%d/%y %H:%M:%S'
echo Moving dilp Read1 File to Read1 Analysis Directory Complete >> ../../../../../analysisnotes/Analysis.log
date '+%m/%d/%y %H:%M:%S' >> ../../../../../analysisnotes/Analysis.log

cd ../read2
#Current Directory=ngs/bwape/inputsequences/hold/dilp/read2
echo Moving dilp Read2 File to Read2 Analysis Directory
date '+%m/%d/%y %H:%M:%S'
echo Moving dilp Read2 File to Read2 Analysis Directory >> ../../../../../analysisnotes/Analysis.log
date '+%m/%d/%y %H:%M:%S' >> ../../../../../analysisnotes/Analysis.log
echo Moving the following files:
for ligne in `ls *.txt`
do                                                                     
echo $ligne
done
for ligne in `ls *.txt`
do
echo Moving the following files: >> ../../../../../analysisnotes/Analysis.log
done
for ligne in `ls *.txt`
do
mv $ligne ../../../illumina/read2/
done
echo Moving dilp Read2 File to Read2 Analysis Directory Complete
date '+%m/%d/%y %H:%M:%S'
echo Moving dilp Read2 File to Read2 Analysis Directory Complete >> ../../../../../analysisnotes/Analysis.log
date '+%m/%d/%y %H:%M:%S' >> ../../../../../analysisnotes/Analysis.log

cd ../../../../../
#Current Directory=ngs/

# Now we call the "BWApe_hg18_v1.sh" script to analyze this sample/lane using bwa sampe

echo starting bwa nanlysis for dilp files using the script /export/Assa/ngs/scripts/BWApe_dm3_V1.sh
/export/Assa/ngs/scripts/BWApe_dm3_V1.sh

echo ***dilp Analysis Complete***

# The analysis directory should now be empty and we can now load the sample2/wt data into the analysis directories

#Check if files exist in the wt hold folder

#List of directoryies to check
temp4=/export/Assa/ngs/bwape/inputsequences/hold/wt/read1
temp5=/export/Assa/ngs/bwape/inputsequences/hold/wt/read2

echo ***Checking wt Hold Folder***
if [ `ls $temp4 | wc -l` != 1 ]       
	then 
	echo " The wt Read1 hold folder does not contain the expect single file "
	echo " ERROR - The script was automatically killed due to a launch error - See Above Error Message"
	exit 2                     
fi
if [ `ls $temp5 | wc -l` != 1 ]       
	then 
	echo " The wt Read2 hold folder does not contain the expect single file "
	echo " ERROR - The script was automatically killed due to a launch error - See Above Error Message"
	exit 2                     
fi
echo ***Found Expected Files***
echo ***Starting The Analysis of wt***
date '+%m/%d/%y %H:%M:%S'
echo ***Starting The Analysis of wt*** >> analysisnotes/Analysis.log
date '+%m/%d/%y %H:%M:%S' >> analysisnotes/Analysis.log

#In the next step we move the "wt" data from ngs/bwape/inputsequences/hold/wt/(read1-2) to ngs/bwape/inputsequences/illumina/(read1-2)

cd bwape/inputsequences/hold/wt/read1
#Current Directory=ngs/bwape/inputsequences/hold/wt/read1
echo Moving wt Read1 File to Read1 Analysis Directory
date '+%m/%d/%y %H:%M:%S'
echo Moving wt Read1 File to Read1 Analysis Directory >> ../../../../../analysisnotes/Analysis.log
date '+%m/%d/%y %H:%M:%S' >> ../../../../../analysisnotes/Analysis.log
echo Moving the following files:
for ligne in `ls *.txt`
do                                                                     
echo $ligne
done
for ligne in `ls *.txt`
do
echo Moving the following files: >> ../../../../../analysisnotes/Analysis.log
done
for ligne in `ls *.txt`
do
mv $ligne ../../../illumina/read1/
done
echo Moving wt Read1 File to Read1 Analysis Directory Complete
date '+%m/%d/%y %H:%M:%S'
echo Moving wt Read1 File to Read1 Analysis Directory Complete >> ../../../../../analysisnotes/Analysis.log
date '+%m/%d/%y %H:%M:%S' >> ../../../../../analysisnotes/Analysis.log

cd ../read2
#Current Directory=ngs/bwape/inputsequences/hold/wt/read2
echo Moving wt Read2 File to Read2 Analysis Directory
date '+%m/%d/%y %H:%M:%S'
echo Moving wt Read2 File to Read2 Analysis Directory >> ../../../../../analysisnotes/Analysis.log
date '+%m/%d/%y %H:%M:%S' >> ../../../../../analysisnotes/Analysis.log
echo Moving the following files:
for ligne in `ls *.txt`
do                                                                     
echo $ligne
done
for ligne in `ls *.txt`
do
echo Moving the following files: >> ../../../../../analysisnotes/Analysis.log
done
for ligne in `ls *.txt`
do
mv $ligne ../../../illumina/read2/
done
echo Moving wt Read2 File to Read2 Analysis Directory Complete
date '+%m/%d/%y %H:%M:%S'
echo Moving wt Read2 File to Read2 Analysis Directory Complete >> ../../../../../analysisnotes/Analysis.log
date '+%m/%d/%y %H:%M:%S' >> ../../../../../analysisnotes/Analysis.log

cd ../../../../../
#Current Directory=ngs/

# Now we call the "BWApe_hg18_v1.sh" script to analyze this sample/lane using bwa sampe

echo starting bwa nanlysis for wt files using the script /export/Assa/ngs/scripts/BWApe_dm3_V1.sh

/export/Assa/ngs/scripts/BWApe_dm3_V1.sh

echo ***wt Analysis Complete***

# The analysis directory should now be empty and we can now load the sample3/lane3 data into the analysis directories

#Check if files exist in the lane3 hold folder

#List of directoryies to check
temp6=/export/Assa/ngs/bwape/inputsequences/hold/lane3/read1
temp7=/export/Assa/ngs/bwape/inputsequences/hold/lane3/read2

echo ***Checking Lane3 Hold Folder***
if [ `ls $temp6 | wc -l` != 1 ]       
	then 
	echo " The Lane3 Read1 hold folder does not contain the expect single file "
	echo " Do you have a third lane? if yes, check the $temp6 directory "
	echo " ERROR - The script was automatically killed due to a launch error - See Above Error Message"
	exit 2                     
fi
if [ `ls $temp7 | wc -l` != 1 ]       
	then 
	echo " The Lane3 Read2 hold folder does not contain the expect single file "
		echo " Do you have a third lane? if yes, check the $temp7 directory "
	echo " ERROR - The script was automatically killed due to a launch error - See Above Error Message"
	exit 2                     
fi
echo ***Found Expected Files***
#Current directory=ngs
echo ***Starting The Analysis of Lane3***
date '+%m/%d/%y %H:%M:%S'
echo ***Starting The Analysis of Lane3*** >> analysisnotes/Analysis.log
date '+%m/%d/%y %H:%M:%S' >> analysisnotes/Analysis.log

#In the next step we move the "Lane3" data from ngs/bwape/inputsequences/hold/lane3/(read1-2) to ngs/bwape/inputsequences/illumina/(read1-2)

cd bwape/inputsequences/hold/lane3/read1
#Current Directory=ngs/bwape/inputsequences/hold/lane3/read1
echo Moving Lane3 Read1 File to Read1 Analysis Directory
date '+%m/%d/%y %H:%M:%S'
echo Moving Lane3 Read1 File to Read1 Analysis Directory >> ../../../../../analysisnotes/Analysis.log
date '+%m/%d/%y %H:%M:%S' >> ../../../../../analysisnotes/Analysis.log
echo Moving the following files:
for ligne in `ls *.txt`
do                                                                     
echo $ligne
done
for ligne in `ls *.txt`
do
echo Moving the following files: >> ../../../../../analysisnotes/Analysis.log
done
for ligne in `ls *.txt`
do
mv $ligne ../../../illumina/read1/
done
echo Moving Lane3 Read1 File to Read1 Analysis Directory Complete
date '+%m/%d/%y %H:%M:%S'
echo Moving Lane3 Read1 File to Read1 Analysis Directory Complete >> ../../../../../analysisnotes/Analysis.log
date '+%m/%d/%y %H:%M:%S' >> ../../../../../analysisnotes/Analysis.log

cd ../read2
#Current Directory=ngs/bwape/inputsequences/hold/lane3/read2
echo Moving Lane3 Read2 File to Read2 Analysis Directory
date '+%m/%d/%y %H:%M:%S'
echo Moving Lane3 Read2 File to Read2 Analysis Directory >> ../../../../../analysisnotes/Analysis.log
date '+%m/%d/%y %H:%M:%S' >> ../../../../../analysisnotes/Analysis.log
echo Moving the following files:
for ligne in `ls *.txt`
do                                                                     
echo $ligne
done
for ligne in `ls *.txt`
do
echo Moving the following files: >> ../../../../../analysisnotes/Analysis.log
done
for ligne in `ls *.txt`
do
mv $ligne ../../../illumina/read2/
done
echo Moving Lane3 Read2 File to Read2 Analysis Directory Complete
date '+%m/%d/%y %H:%M:%S'
echo Moving Lane3 Read2 File to Read2 Analysis Directory Complete >> ../../../../../analysisnotes/Analysis.log
date '+%m/%d/%y %H:%M:%S' >> ../../../../../analysisnotes/Analysis.log

cd ../../../../../
#Current Directory=ngs/

# Now we call the "BWApe_hg18_v1.sh" script to analyze this sample/lane using bwa sampe

/export/Assa/ngs/scripts/BWApe_dm3_V1.sh

echo ***Lane3 Analysis Complete***

# The analysis directory should now be empty and we can now load the sample4/lane4 data into the analysis directories

#Check if files exist in the lane4 hold folder

#List of directoryies to check
temp8=/export/Assa/ngs/bwape/inputsequences/hold/lane4/read1
temp9=/export/Assa/ngs/bwape/inputsequences/hold/lane4/read2

echo ***Checking Lane4 Hold Folder***
if [ `ls $temp8 | wc -l` != 1 ]       
	then 
	echo " The Lane4 Read1 hold folder does not contain the expect single file "
	echo " ERROR - The script was automatically killed due to a launch error - See Above Error Message"
	exit 2                     
fi
if [ `ls $temp9 | wc -l` != 1 ]       
	then 
	echo " The Lane4 Read2 hold folder does not contain the expect single file "
	echo " ERROR - The script was automatically killed due to a launch error - See Above Error Message"
	exit 2                     
fi
echo ***Found Expected Files***
echo ***Starting The Analysis of Lane4***
date '+%m/%d/%y %H:%M:%S'
echo ***Starting The Analysis of Lane4*** >> analysisnotes/Analysis.log
date '+%m/%d/%y %H:%M:%S' >> analysisnotes/Analysis.log

#In the next step we move the "Lane4" data from ngs/bwape/inputsequences/hold/lane4/(read1-2) to ngs/bwape/inputsequences/illumina/(read1-2)

cd bwape/inputsequences/hold/lane4/read1
#Current Directory=ngs/bwape/inputsequences/hold/lane4/read1
echo Moving Lane4 Read1 File to Read1 Analysis Directory
date '+%m/%d/%y %H:%M:%S'
echo Moving Lane4 Read1 File to Read1 Analysis Directory >> ../../../../../analysisnotes/Analysis.log
date '+%m/%d/%y %H:%M:%S' >> ../../../../../analysisnotes/Analysis.log
echo Moving the following files:
for ligne in `ls *.txt`
do                                                                     
echo $ligne
done
for ligne in `ls *.txt`
do
echo Moving the following files: >> ../../../../../analysisnotes/Analysis.log
done
for ligne in `ls *.txt`
do
mv $ligne ../../../illumina/read1/
done
echo Moving Lane4 Read1 File to Read1 Analysis Directory Complete
date '+%m/%d/%y %H:%M:%S'
echo Moving Lane4 Read1 File to Read1 Analysis Directory Complete >> ../../../../../analysisnotes/Analysis.log
date '+%m/%d/%y %H:%M:%S' >> ../../../../../analysisnotes/Analysis.log

cd ../read2
#Current Directory=ngs/bwape/inputsequences/hold/lane4/read2
echo Moving Lane4 Read2 File to Read2 Analysis Directory
date '+%m/%d/%y %H:%M:%S'
echo Moving Lane4 Read2 File to Read2 Analysis Directory >> ../../../../../analysisnotes/Analysis.log
date '+%m/%d/%y %H:%M:%S' >> ../../../../../analysisnotes/Analysis.log
echo Moving the following files:
for ligne in `ls *.txt`
do                                                                     
echo $ligne
done
for ligne in `ls *.txt`
do
echo Moving the following files: >> ../../../../../analysisnotes/Analysis.log
done
for ligne in `ls *.txt`
do
mv $ligne ../../../illumina/read2/
done
echo Moving Lane4 Read2 File to Read2 Analysis Directory Complete
date '+%m/%d/%y %H:%M:%S'
echo Moving Lane4 Read2 File to Read2 Analysis Directory Complete >> ../../../../../analysisnotes/Analysis.log
date '+%m/%d/%y %H:%M:%S' >> ../../../../../analysisnotes/Analysis.log

cd ../../../../../
#Current Directory=ngs/

# Now we call the "BWApe_hg18_v1.sh" script to analyze this sample/lane using bwa sampe

/export/Assa/ngs/scripts/BWApe_dm3_V1.sh

echo ***Lane4 Analysis Complete***

# The analysis directory should now be empty and we can now load the sample5/lane5 data into the analysis directories

#Check if files exist in the dilp hold folder

#List of directoryies to check
temp10=/export/Assa/ngs/bwape/inputsequences/hold/lane5/read1
temp11=/export/Assa/ngs/bwape/inputsequences/hold/lane5/read2

echo ***Checking dilp Hold Folder***
if [ `ls $temp10 | wc -l` != 1 ]       
	then 
	echo " The Lane5 Read1 hold folder does not contain the expect single file "
	echo " ERROR - The script was automatically killed due to a launch error - See Above Error Message"
	exit 2                     
fi
if [ `ls $temp11 | wc -l` != 1 ]       
	then 
	echo " The Lane5 Read2 hold folder does not contain the expect single file "
	echo " ERROR - The script was automatically killed due to a launch error - See Above Error Message"
	exit 2                     
fi
echo ***Found Expected Files***
#Current directory=ngs
echo ***Starting The Analysis of Lane5***
date '+%m/%d/%y %H:%M:%S'
echo ***Starting The Analysis of Lane5*** >> analysisnotes/Analysis.log
date '+%m/%d/%y %H:%M:%S' >> analysisnotes/Analysis.log

#In the next step we move the "Lane5" data from ngs/bwape/inputsequences/hold/lane5/(read1-2) to ngs/bwape/inputsequences/illumina/(read1-2)

cd bwape/inputsequences/hold/lane5/read1
#Current Directory=ngs/bwape/inputsequences/hold/lane5/read1
echo Moving Lane5 Read1 File to Read1 Analysis Directory
date '+%m/%d/%y %H:%M:%S'
echo Moving Lane5 Read1 File to Read1 Analysis Directory >> ../../../../../analysisnotes/Analysis.log
date '+%m/%d/%y %H:%M:%S' >> ../../../../../analysisnotes/Analysis.log
echo Moving the following files:
for ligne in `ls *.txt`
do                                                                     
echo $ligne
done
for ligne in `ls *.txt`
do
echo Moving the following files: >> ../../../../../analysisnotes/Analysis.log
done
for ligne in `ls *.txt`
do
mv $ligne ../../../illumina/read1/
done
echo Moving Lane5 Read1 File to Read1 Analysis Directory Complete
date '+%m/%d/%y %H:%M:%S'
echo Moving Lane5 Read1 File to Read1 Analysis Directory Complete >> ../../../../../analysisnotes/Analysis.log
date '+%m/%d/%y %H:%M:%S' >> ../../../../../analysisnotes/Analysis.log

cd ../read2
#Current Directory=ngs/bwape/inputsequences/hold/lane5/read2
echo Moving Lane5 Read2 File to Read2 Analysis Directory
date '+%m/%d/%y %H:%M:%S'
echo Moving Lane5 Read2 File to Read2 Analysis Directory >> ../../../../../analysisnotes/Analysis.log
date '+%m/%d/%y %H:%M:%S' >> ../../../../../analysisnotes/Analysis.log
echo Moving the following files:
for ligne in `ls *.txt`
do                                                                     
echo $ligne
done
for ligne in `ls *.txt`
do
echo Moving the following files: >> ../../../../../analysisnotes/Analysis.log
done
for ligne in `ls *.txt`
do
mv $ligne ../../../illumina/read2/
done
echo Moving Lane5 Read2 File to Read2 Analysis Directory Complete
date '+%m/%d/%y %H:%M:%S'
echo Moving Lane5 Read2 File to Read2 Analysis Directory Complete >> ../../../../../analysisnotes/Analysis.log
date '+%m/%d/%y %H:%M:%S' >> ../../../../../analysisnotes/Analysis.log

cd ../../../../../
#Current Directory=ngs/

# Now we call the "BWApe_hg18_v1.sh" script to analyze this sample/lane using bwa sampe

/export/Assa/ngs/scripts/BWApe_dm3_V1.sh

echo ***Lane5 Analysis Complete***

# The analysis directory should now be empty and we can now load the sample6/lane6 data into the analysis directories

#Check if files exist in the lane6 hold folder

#List of directoryies to check
temp12=/export/Assa/ngs/bwape/inputsequences/hold/lane6/read1
temp13=/export/Assa/ngs/bwape/inputsequences/hold/lane6/read2

echo ***Checking Lane6 Hold Folder***
if [ `ls $temp12 | wc -l` != 1 ]       
	then 
	echo " The Lane6 Read1 hold folder does not contain the expect single file "
	echo " ERROR - The script was automatically killed due to a launch error - See Above Error Message"
	exit 2                     
fi
if [ `ls $temp13 | wc -l` != 1 ]       
	then 
	echo " The Lane6 Read2 hold folder does not contain the expect single file "
	echo " ERROR - The script was automatically killed due to a launch error - See Above Error Message"
	exit 2                     
fi
echo ***Found Expected Files***
echo ***Starting The Analysis of Lane6***
date '+%m/%d/%y %H:%M:%S'
echo ***Starting The Analysis of Lane6*** >> analysisnotes/Analysis.log
date '+%m/%d/%y %H:%M:%S' >> analysisnotes/Analysis.log

#In the next step we move the "Lane6" data from ngs/bwape/inputsequences/hold/lane6/(read1-2) to ngs/bwape/inputsequences/illumina/(read1-2)

cd bwape/inputsequences/hold/lane6/read1
#Current Directory=ngs/bwape/inputsequences/hold/lane6/read1
echo Moving Lane6 Read1 File to Read1 Analysis Directory
date '+%m/%d/%y %H:%M:%S'
echo Moving Lane6 Read1 File to Read1 Analysis Directory >> ../../../../../analysisnotes/Analysis.log
date '+%m/%d/%y %H:%M:%S' >> ../../../../../analysisnotes/Analysis.log
echo Moving the following files:
for ligne in `ls *.txt`
do                                                                     
echo $ligne
done
for ligne in `ls *.txt`
do
echo Moving the following files: >> ../../../../../analysisnotes/Analysis.log
done
for ligne in `ls *.txt`
do
mv $ligne ../../../illumina/read1/
done
echo Moving Lane6 Read1 File to Read1 Analysis Directory Complete
date '+%m/%d/%y %H:%M:%S'
echo Moving Lane6 Read1 File to Read1 Analysis Directory Complete >> ../../../../../analysisnotes/Analysis.log
date '+%m/%d/%y %H:%M:%S' >> ../../../../../analysisnotes/Analysis.log

cd ../read2
#Current Directory=ngs/bwape/inputsequences/hold/lane6/read2
echo Moving Lane6 Read2 File to Read2 Analysis Directory
date '+%m/%d/%y %H:%M:%S'
echo Moving Lane6 Read2 File to Read2 Analysis Directory >> ../../../../../analysisnotes/Analysis.log
date '+%m/%d/%y %H:%M:%S' >> ../../../../../analysisnotes/Analysis.log
echo Moving the following files:
for ligne in `ls *.txt`
do                                                                     
echo $ligne
done
for ligne in `ls *.txt`
do
echo Moving the following files: >> ../../../../../analysisnotes/Analysis.log
done
for ligne in `ls *.txt`
do
mv $ligne ../../../illumina/read2/
done
echo Moving Lane6 Read2 File to Read2 Analysis Directory Complete
date '+%m/%d/%y %H:%M:%S'
echo Moving Lane6 Read2 File to Read2 Analysis Directory Complete >> ../../../../../analysisnotes/Analysis.log
date '+%m/%d/%y %H:%M:%S' >> ../../../../../analysisnotes/Analysis.log

cd ../../../../../
#Current Directory=ngs/

# Now we call the "BWApe_hg18_v1.sh" script to analyze this sample/lane using bwa sampe

/export/Assa/ngs/scripts/BWApe_dm3_V1.sh

echo ***Lane6 Analysis Complete***

# The analysis directory should now be empty and we can now load the sample7/lane7 data into the analysis directories

#Check if files exist in the lane7 hold folder

#List of directoryies to check
temp14=/export/Assa/ngs/bwape/inputsequences/hold/lane7/read1
temp15=/export/Assa/ngs/bwape/inputsequences/hold/lane7/read2

echo ***Checking Lane7 Hold Folder***
if [ `ls $temp14 | wc -l` != 1 ]       
	then 
	echo " The Lane7 Read1 hold folder does not contain the expect single file "
	echo " ERROR - The script was automatically killed due to a launch error - See Above Error Message"
	exit 2                     
fi
if [ `ls $temp15 | wc -l` != 1 ]       
	then 
	echo " The Lane7 Read2 hold folder does not contain the expect single file "
	echo " ERROR - The script was automatically killed due to a launch error - See Above Error Message"
	exit 2                     
fi
echo ***Found Expected Files***
#Current directory=ngs
echo ***Starting The Analysis of Lane7***
date '+%m/%d/%y %H:%M:%S'
echo ***Starting The Analysis of Lane7*** >> analysisnotes/Analysis.log
date '+%m/%d/%y %H:%M:%S' >> analysisnotes/Analysis.log

#In the next step we move the "Lane7" data from ngs/bwape/inputsequences/hold/lane7/(read1-2) to ngs/bwape/inputsequences/illumina/(read1-2)

cd bwape/inputsequences/hold/lane7/read1
#Current Directory=ngs/bwape/inputsequences/hold/lane7/read1
echo Moving Lane7 Read1 File to Read1 Analysis Directory
date '+%m/%d/%y %H:%M:%S'
echo Moving Lane7 Read1 File to Read1 Analysis Directory >> ../../../../../analysisnotes/Analysis.log
date '+%m/%d/%y %H:%M:%S' >> ../../../../../analysisnotes/Analysis.log
echo Moving the following files:
for ligne in `ls *.txt`
do                                                                     
echo $ligne
done
for ligne in `ls *.txt`
do
echo Moving the following files: >> ../../../../../analysisnotes/Analysis.log
done
for ligne in `ls *.txt`
do
mv $ligne ../../../illumina/read1/
done
echo Moving Lane7 Read1 File to Read1 Analysis Directory Complete
date '+%m/%d/%y %H:%M:%S'
echo Moving Lane7 Read1 File to Read1 Analysis Directory Complete >> ../../../../../analysisnotes/Analysis.log
date '+%m/%d/%y %H:%M:%S' >> ../../../../../analysisnotes/Analysis.log

cd ../read2
#Current Directory=ngs/bwape/inputsequences/hold/lane3/read2
echo Moving Lane7 Read2 File to Read2 Analysis Directory
date '+%m/%d/%y %H:%M:%S'
echo Moving Lane7 Read2 File to Read2 Analysis Directory >> ../../../../../analysisnotes/Analysis.log
date '+%m/%d/%y %H:%M:%S' >> ../../../../../analysisnotes/Analysis.log
echo Moving the following files:
for ligne in `ls *.txt`
do                                                                     
echo $ligne
done
for ligne in `ls *.txt`
do
echo Moving the following files: >> ../../../../../analysisnotes/Analysis.log
done
for ligne in `ls *.txt`
do
mv $ligne ../../../illumina/read2/
done
echo Moving Lane7 Read2 File to Read2 Analysis Directory Complete
date '+%m/%d/%y %H:%M:%S'
echo Moving Lane7 Read2 File to Read2 Analysis Directory Complete >> ../../../../../analysisnotes/Analysis.log
date '+%m/%d/%y %H:%M:%S' >> ../../../../../analysisnotes/Analysis.log

cd ../../../../../
#Current Directory=ngs/

# Now we call the "BWApe_hg18_v1.sh" script to analyze this sample/lane using bwa sampe

/export/Assa/ngs/scripts/BWApe_dm3_V1.sh

echo ***Lane7 Analysis Complete***

# The analysis directory should now be empty and we can now load the sample8/lane8 data into the analysis directories

#Check if files exist in the lane8 hold folder

#List of directoryies to check
temp16=/export/Assa/ngs/bwape/inputsequences/hold/lane8/read1
temp17=/export/Assa/ngs/bwape/inputsequences/hold/lane8/read2

echo ***Checking Lane8 Hold Folder***
if [ `ls $temp16 | wc -l` != 1 ]       
	then 
	echo " The Lane8 Read1 hold folder does not contain the expect single file "
	echo " ERROR - The script was automatically killed due to a launch error - See Above Error Message"
	exit 2                     
fi
if [ `ls $temp17 | wc -l` != 1 ]       
	then 
	echo " The Lane8 Read2 hold folder does not contain the expect single file "
	echo " ERROR - The script was automatically killed due to a launch error - See Above Error Message"
	exit 2                     
fi
echo ***Found Expected Files***
echo ***Starting The Analysis of Lane8***
date '+%m/%d/%y %H:%M:%S'
echo ***Starting The Analysis of Lane8*** >> analysisnotes/Analysis.log
date '+%m/%d/%y %H:%M:%S' >> analysisnotes/Analysis.log

#In the next step we move the "Lane8" data from ngs/bwape/inputsequences/hold/lane8/(read1-2) to ngs/bwape/inputsequences/illumina/(read1-2)

cd bwape/inputsequences/hold/lane8/read1
#Current Directory=ngs/bwape/inputsequences/hold/lane8/read1
echo Moving Lane8 Read1 File to Read1 Analysis Directory
date '+%m/%d/%y %H:%M:%S'
echo Moving Lane8 Read1 File to Read1 Analysis Directory >> ../../../../../analysisnotes/Analysis.log
date '+%m/%d/%y %H:%M:%S' >> ../../../../../analysisnotes/Analysis.log
echo Moving the following files:
for ligne in `ls *.txt`
do                                                                     
echo $ligne
done
for ligne in `ls *.txt`
do
echo Moving the following files: >> ../../../../../analysisnotes/Analysis.log
done
for ligne in `ls *.txt`
do
mv $ligne ../../../illumina/read1/
done
echo Moving Lane8 Read1 File to Read1 Analysis Directory Complete
date '+%m/%d/%y %H:%M:%S'
echo Moving Lane8 Read1 File to Read1 Analysis Directory Complete >> ../../../../../analysisnotes/Analysis.log
date '+%m/%d/%y %H:%M:%S' >> ../../../../../analysisnotes/Analysis.log

cd ../read2
#Current Directory=ngs/bwape/inputsequences/hold/lane8/read2
echo Moving Lane8 Read2 File to Read2 Analysis Directory
date '+%m/%d/%y %H:%M:%S'
echo Moving Lane8 Read2 File to Read2 Analysis Directory >> ../../../../../analysisnotes/Analysis.log
date '+%m/%d/%y %H:%M:%S' >> ../../../../../analysisnotes/Analysis.log
echo Moving the following files:
for ligne in `ls *.txt`
do                                                                     
echo $ligne
done
for ligne in `ls *.txt`
do
echo Moving the following files: >> ../../../../../analysisnotes/Analysis.log
done
for ligne in `ls *.txt`
do
mv $ligne ../../../illumina/read2/
done
echo Moving Lane8 Read2 File to Read2 Analysis Directory Complete
date '+%m/%d/%y %H:%M:%S'
echo Moving Lane8 Read2 File to Read2 Analysis Directory Complete >> ../../../../../analysisnotes/Analysis.log
date '+%m/%d/%y %H:%M:%S' >> ../../../../../analysisnotes/Analysis.log

cd ../../../../../
#Current Directory=ngs/

# Now we call the "BWApe_hg18_v1.sh" script to analyze this sample/lane using bwa sampe

/export/Assa/ngs/scripts/BWApe_dm3_V1.sh

echo ***Lane8 Analysis Complete***
