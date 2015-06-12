source("0_init.R")
do.paired.comparisons.alchemy <- function(expt.cols=NULL,control.cols=NULL,out.dat=NULL,out.tbl=NULL,out.eps=NULL,controls, experimentals) {
  # collect all the data into a big dataframe#
#  list.nodats <- c("1a","2a","3a","4a","5a","6a","7a","8a")#
  num.rep <- length(expt.cols) 
  load(paste(list.nodats[1], ".dat", sep=""))#
  if(is.null(expt.cols) || is.null(control.cols)) {#
    cat("Here are the numbers assigned to your CEL files:\n")#
    print(cbind(1:ncol(expr.norm),colnames(expr.norm)))#
    num.rep <- as.numeric(readline(prompt="How many replicates? "))#
    control.cols <- c()#
    expt.cols <- c()#
    for(i in 1:num.rep) {#
      j <- as.numeric(readline(prompt=paste("Please type the number for control sample in replicate ",i,": ",sep="")))#
      k <- as.numeric(readline(prompt=paste("Please type the number for experimental sample in replicate ",i,": ",sep="")))#
      control.cols <- c(control.cols,j)#
      expt.cols <- c(expt.cols,k)#
    }#
  }#
  if(is.null(out.dat)) {#
     
    out.dat <- readline(prompt="Please type the filename for the result R data file: ")
  }#
  if(is.null(out.tbl)) {
    out.tbl <- readline(prompt="Please type the filename for the result text file: ")
  }#
  if(is.null(out.eps)) {#
    out.eps <- readline(prompt="Please type the filename for the result image file: ")
  }#
#
 
  cat("Performing a comparison for paired samples:\n")
  tmp.show <- cbind(Expt=colnames(expr.norm)[expt.cols],Control=colnames(expr.norm)[control.cols])#
  print(tmp.show)#
  cat("Writing results to ",out.dat,", ",out.tbl,", and ",out.eps,".\n",sep="")
#
  num.pseudorep <- length(list.nodats)
 
  tmp.conf <- floor(10 * (num.pseudorep * num.rep - 1) / (num.rep - 1))
  controls <- CONTROL.ARRAYS
  experimentals <- TREATMENT.ARRAYS

  res <- paired.composite(list.nodats,control.cols,expt.cols,conf=tmp.conf)#
  save(res,file=out.dat)
  write.table(cbind(res$d.orig,res$mean.sig,res$all.signal),file=out.tbl,sep="\t",quote=FALSE)
  postscript(file=out.eps,horizontal=FALSE,onefile=FALSE,height=5,width=5,points=12)
  par(mar=c(4,4,0.5,0.5))
  q.cutoffs <- c(0.00001,0.0001,0.001,0.002,0.005,0.01,0.02,0.05,0.1)
  num.positive <- c()
  for (i in q.cutoffs) {
    num.positive <- c(num.positive,sum(res$d.orig$q <= i))
  }
  plot(log10(q.cutoffs),num.positive,pch=19,xlab="Log(10) Q-value cutoff",ylab="Number called differentially expressed")
  lines(log10(q.cutoffs),num.positive)
  dev.off()
}


library(affy)
library(goldenspike)
library(gcrma)

list.nodats <- c("2_normalized/2nd/mas5_loess_loess","2_normalized/2nd/mas5_loess_vsn","2_normalized/2nd/mas5_quantiles_loess","2_normalized/2nd/mas5_quantiles_vsn",
"2_normalized/2nd/gc_loess_loess","2_normalized/2nd/gc_loess_vsn","2_normalized/2nd/gc_quantiles_loess","2_normalized/2nd/gc_quantiles_vsn")

do.paired.comparisons.alchemy(expt.cols=TREATMENT.ARRAYS,control.cols=CONTROL.ARRAYS, out.dat="3_simplePierreAlchemy/paired_result_present.dat",
out.tbl="3_simplePierreAlchemy/paired_result_present.tbl", out.eps="3_simplePierreAlchemy/paired_result_present.eps")

