
make.expr.summaries.alchemy <- function(filenames=character(0)) {
 
  require(affy)
  require(gcrma)
  require(vsn)

  save(data,file="data.orig")
  option<-commandArgs()[length(commandArgs())]


   if (option!="is2ndOnly"){
     cat("Run first normalization");cat("\n")

# calculate MAS5 epression values
  load("data.orig")
  expression <- expresso(data,bg.correct=TRUE,bgcorrect.method="mas",
      normalize=TRUE,normalize.method="loess", normalize.param = list(span=1/10,subset = sample(1:(dim(pm(data))[1]), min(c(50000, nrow(pm(data))))) ),
      pmcorrect.method="mas",
      summary.method="medianpolish",
      verbose=TRUE)
  print("finised loess normalization with mas5 pm correction. Saving to 1st/mas5_loess.dat")
  print(paste("thw working directory is: ", getwd()))
  write.table(exprs(expression), file="1st/exprsVal_mas5_loess.txt", sep="\t", quote=F)
  save(expression,file="1st/mas5_loess.dat")
  
  load("data.orig")
  expression <- expresso(data,bg.correct=TRUE,bgcorrect.method="mas",
      normalize=TRUE,normalize.method="quantiles",
      pmcorrect.method="mas",
      summary.method="medianpolish",
      verbose=TRUE)
  print("finised quantile normalization with mas5 pm correction. Saving to 1st/mas5_quantiles.dat")
  print(paste("thw working directory is: ", getwd()))
  write.table(exprs(expression), file="1st/exprsVal_mas5_quantile.txt", sep="\t", quote=F)
  save(expression,file="1st/mas5_quantiles.dat")
  


    #   GC-RMA transform the PM probe values using the MM probes
    #   as the negative controls for NSB correction
    load("data.orig")
    data.gcrma <- bg.adjust.gcrma(data, affinity.source="local", type="fullmodel", fast=FALSE)
    save(data.gcrma, file = "data_gcrma.dat")
      

  load("data_gcrma.dat")
 #keep present calls only for second round; the ﬁrst argument subsets the features and the second argument subsets the samples.
    expression <- expresso(data.gcrma,bg.correct=FALSE,
      normalize=TRUE,normalize.method="loess", normalize.param = list(span=1/10,subset = sample(1:(dim(pm(data))[1]), min(c(50000, nrow(pm(data))))) ),
      pmcorrect.method="pmonly",
      summary.method="medianpolish",
      verbose=TRUE)
  print("finised loess normalization with pm_only (gcrma) pm correction. Saving to 1st/gc_loess.dat")
  write.table(exprs(expression), file="1st/exprsVal_gc_loess.txt", sep="\t", quote=F)
  save(expression,file="1st/gc_loess.dat")
  
  load("data_gcrma.dat")
  expression <- expresso(data.gcrma,bg.correct=FALSE,
      normalize=TRUE,normalize.method="quantiles",
      pmcorrect.method="pmonly",
      summary.method="medianpolish",
      verbose=TRUE)
  print("finised quantile normalizationwith pm_only (gcrma) pm correction. Saving to 1st/gc_quantiles.dat")
  write.table(exprs(expression), file="1st/exprsVal_gc_quantile.txt", sep="\t", quote=F)
  save(expression,file="1st/gc_quantiles.dat")
   }
  cat("finished 1st normalization"); cat("\n"); cat("\n"); cat("\n")
  ######################################################
  ## now perform a second loess on these datasets.    ##
  ######################################################
  list.nodats <- c("mas5_loess","mas5_quantiles",
  "gc_loess","gc_quantiles")

          cat("Run 2nd normalization");cat("\n");cat("\n");cat("\n")

    for(i in list.nodats) {
        filename <- paste("1st/", i,".dat",sep="")
        load(filename)
	    #remove absent probesets before the second round of normalization!!!
        # you can't remove them for first round: needed for non-specific back ground correction
	    # #################################

        allfeatures<-length(featureNames(expression))
	    expression<-expression[as.vector(present),]
        presentfeatures<-length(featureNames(expression))
	    if (exists("GENE.FILTER.2ND.NORM")) {
	        print(paste(length(GENE.FILTER.2ND.NORM), "are known to be in the ovary. These genes should be excluded"));cat("\n")
	        tmp<-GENE.FILTER.2ND.NORM[GENE.FILTER.2ND.NORM%in%present]
                print(paste(length(tmp), "genes are were found present and will be removed"));cat("\n")	        
                expression<-expression[!featureNames(expression)%in% as.vector(tmp),]
	    }

        remainingfeatures<-length(featureNames(expression))

        cat("All features: ") ;cat(allfeatures) ; cat("\n");
        cat("Present features: ") ;cat(presentfeatures) ; cat("\n");
        cat("Remaining features after filtering: ") ;cat(remainingfeatures) ;cat("\n");

        expr <- prepare.expr(exprs(expression),log.it=TRUE,cutoff.value=-3)

        gene.names <- rownames(expr)
        expr.norm <- normalize.loesssubset(expr, subset = 1:length(gene.names),log.it=FALSE,span=1/10,sample.length=nrow(expr))
        expr.norm <- 2^expr.norm
        s <- 500/mean(expr.norm, trim=0.02)
        expr.norm <- s*expr.norm
        expr.norm <- prepare.expr(expr.norm,log.it=TRUE,cutoff.value=-3)
        filename <- paste("2nd/",i,"_loess.dat",sep="")
        write.table(expr.norm, file=paste(filename, "_loess.txt", sep=""), sep="\t", quote=F)
        save(expr.norm, file=filename)
        rm(expr.norm)
        filename <- paste("1st/",i,".dat",sep="")
        load(filename)
	    #remove absent probesets before the second round of normalization!!!
        # you can't remove them for first round: needed for non-specific back ground correction
	    # #################################
        expression<-expression[as.vector(present),]
        if (exists("GENE.FILTER.2ND.NORM")) {
	        print(paste(length(GENE.FILTER.2ND.NORM), "are known to be in the ovary. These genes should be excluded"));cat("\n")
	        tmp<-GENE.FILTER.2ND.NORM[GENE.FILTER.2ND.NORM%in%present]
	        print(paste(length(tmp), "genes are were found present and will be removed"));cat("\n")
	        expression<-expression[!featureNames(expression)%in% as.vector(tmp),]
	    }
        expr <- prepare.expr(exprs(expression),log.it=TRUE,cutoff.value=-3)
        expr <- 2^expr
        vsn2(expr, lts.quantile = 0.85) -> expr.norm
        exprs(expr.norm) -> expr.norm
        expr.norm <- exp(expr.norm)
        s <- 500/mean(expr.norm, trim=0.02)
        expr.norm <- s*expr.norm
        expr.norm <- prepare.expr(expr.norm,log.it=TRUE,cutoff.value=-3)
        filename <- paste("2nd/",i,"_vsn.dat",sep="")
        write.table(expr.norm, file=paste(filename, "_vsn.txt", sep=""), sep="\t", quote=F)
        save(expr.norm, file=filename)
        rm(expr.norm)
  }
}







# Same as in 4a except that the absent probesets/transcripts are being removed before the second round
# of normalization is done

library(goldenspike)

source("0_init.R")
callsdir = paste(getwd(), "1_calls/affy", sep="/")
rawdir = paste(getwd(), "0_raw", sep="/")
normadir=paste(getwd(), "2_normalized/", sep="/")
print(paste("the working directory now is: ",getwd()))
setwd(callsdir)
present <-read.delim("present_gcrma.txt", header =FALSE, sep="\t")[,1]
print(paste("the number of present IDs is ", length(present)));cat("\n")# number pf genes passing the present-call test
print(paste("the working directory now is: ",getwd()))
setwd(rawdir)
files<-list.files()
celfiles<-files[grep(".CEL",files)]
print(paste("the working directory now is: ",getwd()))
setwd(normadir)
print(paste("the working directory now is: ",getwd()))
make.expr.summaries.alchemy(filenames=celfiles)



