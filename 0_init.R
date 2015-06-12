##!## Todo: move variables that are used in one script only to the corresponding script!

###########################################################################################
####### load mandatory variables
###########################################################################################
BASEDIR<-Sys.getenv('BASEDIR')
SVN.PATH<-Sys.getenv('SVN_PATH')
CONTROL.ARRAYS.FILENAMES<-as.vector(strsplit(Sys.getenv('CONTROL_ARRAYS'), ";")[[1]])
TREATMENT.ARRAYS.FILENAMES<-as.vector(strsplit(Sys.getenv('TREATMENT_ARRAYS'), ";")[[1]])

cat(length(TREATMENT.ARRAYS.FILENAMES),"\n")
cat(TREATMENT.ARRAYS.FILENAMES,"\n")
cat(length(CONTROL.ARRAYS.FILENAMES),"\n")
cat(CONTROL.ARRAYS.FILENAMES,"\n")

EASE.CONFIG.FILE<-Sys.getenv('EASE_CONFIG_FILE')
MODEL_library<- as.character(Sys.getenv('MODEL_library'))

###########################################################################################
####### Load data
###########################################################################################
cat("Load cel files","\n")
# changed to fit the current directory structure.
#setwd(paste(BASEDIR,"0_raw", sep="/"))
setwd("0_raw")
library(affy)
data<-ReadAffy(filenames=c(as.character(TREATMENT.ARRAYS.FILENAMES), as.character(CONTROL.ARRAYS.FILENAMES)))
cat("Number arrays loaded: ", length(data), "\n")
setwd(BASEDIR)

###########################################################################################
####### load additional variables or take default values
###########################################################################################
tmp<-Sys.getenv('GENE_FILTER_2ND_NORM')

if (nchar(Sys.getenv('GENE_FILTER_2ND_NORM'))<1){
  GENE.FILTER.2ND.NORM<- vector()
}else{
cat("Sys.getenv('GENE_FILTER_2ND_NORM')");
  GENE.FILTER.2ND.NORM<-as.vector(read.table(Sys.getenv('GENE_FILTER_2ND_NORM'), sep="\n")[,1])
}

CONTROL.ARRAYS<-which(sampleNames(data) %in% CONTROL.ARRAYS.FILENAMES)
TREATMENT.ARRAYS<-which(sampleNames(data) %in% TREATMENT.ARRAYS.FILENAMES)

cat("Control Arrays: ", CONTROL.ARRAYS, "\n")
cat("Control Arrays: ",TREATMENT.ARRAYS, "\n")

MINIMUM<-as.double(Sys.getenv('MINIMUM'))
MAXIMUM<-as.double(Sys.getenv('MAXIMUM'))
STEPS<-as.double(Sys.getenv('STEPS'))

###########################################################################################
####### Print settings to log file
###########################################################################################


cat("#####################################################################")
cat("\n")
cat("##############      Settings             ############################")
cat("\n")
cat("#####################################################################")
cat("\n")
cat("BASEDIR:");cat("\t");cat(BASEDIR);cat("\n")
cat("SVN.PATH:");cat("\t");cat(SVN.PATH);cat("\n")
cat("CONTROL.ARRAYS:");cat("\t");cat(CONTROL.ARRAYS);cat("\n")
cat("TREATMENT.ARRAYS:");cat("\t");cat(TREATMENT.ARRAYS);cat("\n")
cat("MODEL_library:");cat("\t");cat(MODEL_library);cat("\n")
cat("MINIMUM:");cat("\t");cat(MINIMUM);cat("\n")
cat("MAXIMUM:");cat("\t");cat(MAXIMUM);cat("\n")
cat("STEPS:");cat("\t");cat(STEPS);cat("\n")
cat("EASE.CONFIG.FILE:");cat("\t");cat(EASE.CONFIG.FILE);cat("\n")
cat("GENE.FILTER.2ND.NORM");cat("\n");cat(length(GENE.FILTER.2ND.NORM));cat("\n")

print("#####################################################################")
print("\n")

print("#####################################################################")
print("\n")
print("##############      Settings             ############################")
print("\n")
print("#####################################################################")
print("\n")
print("BASEDIR:");print("\t");print(BASEDIR);print("\n")
print("SVN.PATH:");print("\t");print(SVN.PATH);print("\n")
print("CONTROL.ARRAYS:");print("\t");print(CONTROL.ARRAYS);print("\n")
print("TREATMENT.ARRAYS:");print("\t");print(TREATMENT.ARRAYS);print("\n")
print("MODEL_library:");print("\t");print(MODEL_library);print("\n")
print("MINIMUM:");print("\t");print(MINIMUM);print("\n")
print("MAXIMUM:");print("\t");print(MAXIMUM);print("\n")
print("STEPS:");print("\t");print(STEPS);print("\n")
print("EASE.CONFIG.FILE:");print("\t");print(EASE.CONFIG.FILE);print("\n")
print("GENE.FILTER.2ND.NORM");print("\n");print(length(GENE.FILTER.2ND.NORM));print("\n")

print("#####################################################################")
print("\n")

