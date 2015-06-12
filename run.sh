#/bin/bash

##########################################
# Required Command Line Input Params     #
##########################################
MINPARAMS=1
CONFIG_FILE=$1

##########################################
# Load Variables given in the Config file#
##########################################
source $CONFIG_FILE
echo "`tput bold`the config file is at:`tput sgr0` $CONFIG_FILE"
echo "`tput bold`my working directory is at:`tput sgr0` $BASEDIR"
if [ $# -lt "$MINPARAMS" ]
then
  echo
  echo "This script needs at least $MINPARAMS command-line arguments!"
fi  

if [ -z "$CONFIG_FILE" ]
then
		echo "Enter path for the configuration file: "
		read CONFIG_FILE
		
		if test ! -r "$CONFIG_FILE" ; then
		    echo "$CONFIG_FILE does not exist or is not readable."
		    exit 1
		fi
fi

echo "`tput bold`Are the config_pipeline paths changed?`tput sgr0`"
echo "`tput bold`Change the path to the working directory and the names of the CEL files!!!`tput sgr0`"

echo
echo "`tput bold`Run full analysis?`tput sgr0` `tput bold`(y)`tput sgr0`es, `tput bold`(n)`tput sgr0`o: "
echo "to run CatMap, please choose option `tput bold`e`tput sgr0` first to preper the category lists  "



read -e RUNPARAM

if [ "$RUNPARAM" == "n" ]; then
	echo
	echo "`tput bold`Select STEPS to run:`tput sgr0`"
	echo "
	mkdir(1),
	plots(2),
	calls(3),
	normalization(4),
	normalization_2ndonly(f),	
	simplePierre(5),
	exports(6),
	ease(7),
	getMartSequences(8),
	clover(9),
	clover_post(a),
	catmap(b),
	catmap_prepare(e)
	limma(l)
	"
	read -e STEPPARAM

else
	echo "Running all STEPS ..."
	# Note that "a" and "d" and "f" are exlcuded because these are postprocessing STEPS that require the data produced by previous STEPS
	STEPPARAM="1234567891011bceghijklmnopqrstuvwxyz"
fi


logDir=$BASEDIR"/10_Logging"


if echo $STEPPARAM | grep "1"
then
   START=$(date +%s)

#######################################################################
##############    Prepare Directory structures           ############## 
#######################################################################
mkdir -p $BASEDIR"/1_calls"
mkdir -p $BASEDIR"/1_calls/affy"
mkdir -p $BASEDIR"/1_calls/geneids"
mkdir -p $BASEDIR"/1_calls/martSequences"
mkdir -p $BASEDIR"/2_normalized"
mkdir -p $BASEDIR"/2_normalized/1st"
mkdir -p $BASEDIR"/2_normalized/2nd"
mkdir -p $BASEDIR"/3_simplePierreAlchemy"
mkdir -p $BASEDIR"/3_simplePierreAlchemy/fractions"
mkdir -p $BASEDIR"/4_genelists"
mkdir -p $BASEDIR"/4_genelists/affy"
mkdir -p $BASEDIR"/4_genelists/geneids"
mkdir -p $BASEDIR"/4_genelists/martSequences"
mkdir -p $BASEDIR"/5_EASE"
mkdir -p $BASEDIR"/5_EASE/full"
mkdir -p $BASEDIR"/5_EASE/slim"
mkdir -p $BASEDIR"/6_Clover"
mkdir -p $BASEDIR"/6_Clover/details"
mkdir -p $BASEDIR"/6_Clover/summary"
mkdir -p $BASEDIR"/7_CatMap"
mkdir -p $BASEDIR"/7_CatMap/result"
mkdir -p $BASEDIR"/10_Logging"
mkdir -p $BASEDIR"/11_lm"

    END=$(date +%s)
    DIFF=$(( $END - $START ))
    echo "It took $DIFF seconds to prepare the directory structure" >> $logDir"/benchmark.log"

fi

set -o errexit  # Exit when one of the following commands fails.  Same as `set -e'


#######################################################################
##############    CALL R                                 ############## 
#######################################################################



if echo $STEPPARAM | grep "2"
then
    START=$(date +%s)
    # run R
    PLOTPARAM=FALSE
    export PLOTPARAM
    R --slave < scripts/plots.R > $logDir"/plots.log"
    END=$(date +%s)
    DIFF=$(( $END - $START ))
    echo "It took $DIFF seconds to produce the raw data plots " >> $logDir"/benchmark.log"
fi

if echo $STEPPARAM | grep "3"
then

    PLOTPARAM=TRUE
    export PLOTPARAM
    START=$(date +%s)
    # run R
    R --slave < scripts/calls.R > $logDir"/calls.log"
    END=$(date +%s)
    DIFF=$(( $END - $START ))
    echo "It took $DIFF seconds to do present/absent calls " >> $logDir"/benchmark.log"
fi

if echo $STEPPARAM | grep "4"
then
    START=$(date +%s)
    # run R
    R --slave --args isBoth < scripts/normalize.R > $logDir"/normalize.log"
    END=$(date +%s)
    DIFF=$(( $END - $START ))
    echo "It took $DIFF seconds to run normalization " >> $logDir"/benchmark.log"
fi

if echo $STEPPARAM | grep "f"
then
    START=$(date +%s)
    # run R
    R --slave --args is2ndOnly < scripts/normalize.R > $logDir"/normalize.log"
    END=$(date +%s)
    DIFF=$(( $END - $START ))
    echo "It took $DIFF seconds to run second step of normalization " >> $logDir"/benchmark.log"
fi



if echo $STEPPARAM | grep "5"
then
    START=$(date +%s)
    # run R

    if ($IS_PAIRED_ANALYSIS)
    then
         R --slave < scripts/simplePierre_paired.R > $logDir"/simplePierre.log"
    else
         R --slave < scripts/simplePierre.R > $logDir"/simplePierre.log"
    fi

    END=$(date +%s)
    DIFF=$(( $END - $START ))
    echo "It took $DIFF seconds to run simplePierre" >> $logDir"/benchmark.log"

fi

if echo $STEPPARAM | grep "6"
then
    START=$(date +%s)
    # run R
    R --slave < scripts/exports.R > $logDir"/exports.log"
    END=$(date +%s)
    DIFF=$(( $END - $START ))
    echo "It took $DIFF seconds to export genelists " >> $logDir"/benchmark.log"

#echo $?
fi



if echo $STEPPARAM | grep "7"
then
    START=$(date +%s)
    # run R
    R --slave < scripts/EASE.R> $logDir"/ease.log"
    END=$(date +%s)
    DIFF=$(( $END - $START ))
    echo "It took $DIFF seconds to run EASE  on genelists" >> $logDir"/benchmark.log"
fi

if echo $STEPPARAM | grep "l"
then
    START=$(date +%s)
    # run R
    R --slave < scripts/limma.R> $logDir"/limma.log"
    END=$(date +%s)
    DIFF=$(( $END - $START ))
    echo "It took $DIFF seconds to run LIMMA  on genelists" >> $logDir"/benchmark.log"
fi


if echo $STEPPARAM | grep "8"
then
    START=$(date +%s)
    # run R
    # get Clover input sequences
    R --slave < scripts/getMartSequences.R > $logDir"/getMartSequences.log"
    END=$(date +%s)
    DIFF=$(( $END - $START ))
    echo "It took $DIFF seconds to get promoter sequences for genelists " >> $logDir"/benchmark.log"

fi


# Requires Research Queue
if echo $STEPPARAM | grep "9"
then
	sh $SVN_PATH"/Pipeline/clover.sh" $BASEDIR $SVN_PATH
    echo "It took ? seconds to run Clover " >> $logDir"/benchmark.log"

fi

# Postprocessing
if echo $STEPPARAM | grep "a"
then
    START=$(date +%s)
    # run R
	echo $SVN_PATH
	echo $BASEDIR
	echo $TF_FACTOR
    java -classpath $SVN_PATH"/Java/classes/" main.util.CloverParserPremium $BASEDIR"/6_Clover/details/output" $BASEDIR"/6_Clover" $TF_FACTOR
    END=$(date +%s)
    DIFF=$(( $END - $START ))
    echo "It took $DIFF seconds to run Clover postprocessing " >> $logDir"/benchmark.log"
fi

#########################################################
######## Still testing mode (contains hard coded paths)
#########################################################
if echo $STEPPARAM | grep "e"
then
    START=$(date +%s)
    # run R
    # get Clover input sequences
    R --slave < scripts/getPresentCategories.R > $logDir"/getPresentCategories.log"
    END=$(date +%s)
    DIFF=$(( $END - $START ))
    echo "It took $DIFF seconds to remove categories that are not present" >> $logDir"/benchmark.log"

fi


if echo $STEPPARAM | grep "b"
then

    	START=$(date +%s)
	file=$BASEDIR"/4_genelists/geneids/unpaired_listfile_increasing.txt"
	if [ -e $file ] ;
        then
	    perl $SVN_PATH"/3rd_party_software/Catmap/Catmap.pl" --categoryfile $BASEDIR"/7_CatMap/category2genes_present.txt" --listfile $BASEDIR"/4_genelists/geneids/unpaired_listfile_increasing.tsv" --randomnull --outputfile $BASEDIR"/7_CatMap/result/result_unpaired_increasing.txt" --companionfile $BASEDIR"/7_CatMap/result/companion_unpaired_increasing.txt"
        echo $?
    fi

    file="$BASEDIR/4_genelists/geneids/unpaired_listfile_decreasing.txt"
	if [ -e $file ] ;
        then
        perl $SVN_PATH"/3rd_party_software/Catmap/Catmap.pl" --categoryfile $BASEDIR"/7_CatMap/category2genes_present.txt" --listfile $BASEDIR"/4_genelists/geneids/unpaired_listfile_decreasing.tsv" --randomnull --outputfile $BASEDIR"/7_CatMap/result/result_unpaired_decreasing.txt" --companionfile $BASEDIR"/7_CatMap/result/companion_unpaired_decreasing.txt"
    fi

    END=$(date +%s)    
    DIFF=$(( $END - $START ))
    echo "It took $DIFF seconds to run Catmap " >> $logDir"/benchmark.log"
fi







