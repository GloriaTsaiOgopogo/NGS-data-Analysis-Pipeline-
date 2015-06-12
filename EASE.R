################################################################################################################# 
################## Functional analysis ##########################################################################
# Reimplemtation of the DAVID/EASE software (http://david.abcc.ncifcrf.gov/ease/ease1.htm)
# 1. advantages: programmatic access, newer annotation files than provided in DAVID/EASE
# 2. disadvantages: annotation files need to be prepared and updated regularly 
#################################################################################################################  
#################################################################################################################  


source("0_init.R")
library(stringr)
setwd(paste(getwd(),"/5_EASE", sep=""))
source(paste(SVN.PATH, "/R/ease/ease_code.R", sep=""))

####################################################################################################
############                               All lists                                    ############ 
####################################################################################################

# Call the method: runEASE(geneList, background, config)
genelists<-list.files("../4_genelists/geneids/")
genelists<-genelists[grep(".txt", genelists)]

for (i in 1:length(genelists)){
cat(genelists[i])	;cat("\n")
if (file.info(paste('../4_genelists/geneids/', genelists[i], sep=''))[1, 1] != 0) {
	
	if (length(grep("regulated", genelists[i]))==0){
		next;
	}
	result <- runEASE(geneList.affyids.file=paste("../4_genelists/geneids",genelists[i],sep="/"), bg.affyids.file="../1_calls/geneids/present_gcrma_geneids.txt", config.file=EASE.CONFIG.FILE)

	# Save the result
	write.table(result, file=paste("full",genelists[i], sep="/"), sep="\t", quote=F, col.names=T, row.names=F) 

	result_slim<-result[result$EASE<0.05,]
	result_slim$category<-substring(as.character(result_slim$category), 0,50)
	result_slim<-result_slim[,c(2,7,9)]


    currentResult<-result_slim

    # rounds the results and reducesss numbers after point.
    for (k in 1:length(currentResult[,2])){
	    for (j in 2:3){
	if (is.na(currentResult[k,j]==TRUE)){
		print(paste("NA detected in line ", k ," at position ", j, " breaking and going to next one", sep=""))
		cat('\n\n\n\n')
		break;
		}
	    a<-currentResult[k,j]

        if (length(grep("e",a))>0){
		    # get exponent
		    begin=str_locate_all(a, "e-")[[1]][2]
		    end=nchar(a)
		    exponent=substring(a, begin+1, end)
		    a=a*10^as.numeric(exponent)
		    a=as.numeric(as.character(round(a, 2)))
		    a=a/10^as.numeric(exponent)
	    }else{
	    	a=as.numeric(as.character(round(a,4)))
	    }
		currentResult[k,j]<-a
		}
	}
    result_slim =currentResult
	write.table(result_slim, file=paste("slim",genelists[i], sep="/"),, sep="\t", quote=F, col.names=T, row.names=F)    
	}
else (cat(paste('going to next one', '\n\n\n', sep='')))
}


