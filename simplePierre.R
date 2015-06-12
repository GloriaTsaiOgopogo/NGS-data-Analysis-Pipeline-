source("0_init.R")
source("updated_simple_pierre.R") # the pierre.unpair.combine methods were updated
print("I run the unpaired analysis of simeplePierre \n")
do.unpaired.comparisons.alchemy <- function(expt.cols=NULL,control.cols=NULL,out.dat=NULL, mean.tbl=NULL, out.tbl=NULL,out.eps=NULL) {

  # collect all the data into a big dataframe
  load(paste(list.nodats[1], ".dat", sep="") )
  if(is.null(expt.cols) || is.null(control.cols)) {
    cat("Here are the numbers assigned to your CEL files:\n")
    print(cbind(1:ncol(expr.norm),colnames(expr.norm)))
    control.cols <- as.numeric(unlist(strsplit(readline(prompt="Please type the numbers for the control samples, separated by spaces: "),'[[:space:]]+')))
    expt.cols <- as.numeric(unlist(strsplit(readline(prompt="Please type the numbers for the experimental samples, separated by spaces: "),'[[:space:]]+')))
  }
  if(is.null(out.dat)) {
  cat(out.dat);cat("\n"); cat(is.null(out.dat))
	out.dat <- readline(prompt="Please type the filename for the result R data file: (without \"\") ")
  }
  if(is.null(out.tbl)) {
    out.tbl <- readline(prompt="Please type the filename for the result text file: (without \"\") ")
  }
  if(is.null(out.eps)) {
    out.eps <- readline(prompt="Please type the filename for the result image file: (without \"\") ")
  }
  
  cat("Performing a comparison for unpaired samples:\n")
  cat("Experimental sample names:\n")
  print( colnames(expr.norm[,expt.cols]))
  cat("\nControl sample names:\n")
  print(colnames(expr.norm[,control.cols]))
  cat("Writing results to ",out.dat,".1, and ",out.dat,".2, \n writing the table of genes to ",out.tbl,"\n  and  the image to the file:",out.eps,".\n",sep="")

  all.expr <- c()
  gene.names <- rownames(expr.norm)
  order.template <- order(gene.names)
  for(i in 1:length(list.nodats)) {
    filename <- paste(list.nodats[i],".dat",sep="")
    load(filename)
    expr <- prepare.expr(expr.norm,log.it=TRUE,cutoff.value=-3)
    rm(expr.norm)
    gene.names <- rownames(expr)
    if(length(which(order(gene.names) != order.template)) > 0) {
      cat("Error: gene names are not in same order for the individual datasets\n")
    }
    all.expr <- c(all.expr,list(i=expr))
  }

  tmp.conf <- 2 * (length(control.cols) + length(expt.cols))
  res <- pierre.unpair.combine.2(all.expr,controls=control.cols,expts=expt.cols,winsize=100,conf=tmp.conf,minrep=3,num.permutations=100)
  save(res,file=paste(out.dat,".1",sep=""))

  # get the mean signal level
  all.sig <- c()
  for(i in 1:length(list.nodats)) {
    tmp.sig <- apply(all.expr[[i]],1,mean.na)
    all.sig <- cbind(all.sig,tmp.sig)
  }
  all.sig <- apply(all.sig,1,mean.na)

  res.summary <- cbind(res$d.orig$combined.mean,res$d.orig$var.a,res$d.orig$var.b,res$d.orig$stat,res$d.orig$q,all.sig, res$d.orig$mean.ctrl, res$d.orig$mean.expt)
  colnames(res.summary) <- c("mean.fc","var.a","var.b","stat","q","mean.sig", "mean.ctrl", "mean.expt")
	#updated part
	#the mean of the two conditions for this fraction will be saved in the 3_simplePierreAlchemy/fractions/ directory
	write.table(as.data.frame(res.summary[,7]), file=paste(mean.tbl, "mean.ctrl", sep="/"), quote=F)
	write.table(as.data.frame(res.summary[,8]), file=paste(mean.tbl, "mean.expt", sep="/"), quote=F)

  save(res.summary,file=paste(out.dat,".2",sep=""))
  write.table(res.summary,file=out.tbl,quote=FALSE,sep="\t")

  postscript(file=out.eps,horizontal=FALSE,onefile=FALSE,height=5,width=5,points=12)
  par(mar=c(4,4,0.5,0.5))
  q.cutoffs <- c(0.00001,0.0001,0.001,0.002,0.005,0.01,0.02,0.05,0.1, 0.5, 0.7, 1)
  num.positive <- c()
  for (i in q.cutoffs) {
    num.positive <- c(num.positive,sum(res.summary[,"q"] <= i))
  }
  plot(log10(q.cutoffs),num.positive,pch=19,xlab="Log(10) Q-value cutoff",ylab="Number called differentially expressed")
  lines(log10(q.cutoffs),num.positive)
  dev.off()
}


library(affy)
library(goldenspike)
library(gcrma)

#list.nodats<-c("2_normalized/rma/rma")
list.nodats <- c("2_normalized/2nd/mas5_loess_loess","2_normalized/2nd/mas5_loess_vsn","2_normalized/2nd/mas5_quantiles_loess","2_normalized/2nd/mas5_quantiles_vsn",
"2_normalized/2nd/gc_loess_loess","2_normalized/2nd/gc_loess_vsn","2_normalized/2nd/gc_quantiles_loess","2_normalized/2nd/gc_quantiles_vsn")

do.unpaired.comparisons.alchemy(expt.cols=TREATMENT.ARRAYS,control.cols=CONTROL.ARRAYS, out.dat="3_simplePierreAlchemy/unpaired_result_present.dat", mean.tbl="3_simplePierreAlchemy/fractions/", 
out.tbl="3_simplePierreAlchemy/unpaired_result_present.tbl", out.eps="3_simplePierreAlchemy/unpaired_result_present.eps")
#expt.cols=TREATMENT.ARRAYS;control.cols=CONTROL.ARRAYS; out.dat="3_simplePierreAlchemy/unpaired_result_present.dat"; mean.tbl="3_simplePierreAlchemy/fractions/"; out.tbl="3_simplePierreAlchemy/unpaired_result_present.tbl"; out.eps="3_simplePierreAlchemy/unpaired_result_present.eps"


