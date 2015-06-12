 
##################################################################################################
###### present_mas5: implementation by Affymetrix
###### 
##################################################################################################
present_mas5<-function(data){
	
	raw.calls <- mas5calls(data,alpha1=.04, alpha2=.06)	
	CONTROL.ARRAYS<-assayData(raw.calls)[["se.exprs"]][,CONTROL.ARRAYS]
	mutantArrays<-assayData(raw.calls)[["se.exprs"]][,TREATMENT.ARRAYS]
	
	control.calls <- apply(CONTROL.ARRAYS, 1, function(x) median(x))
	het.calls <- apply(mutantArrays, 1, function(x) median(x))
	
	present <- all.ps[control.calls < .06  | het.calls  < .06 ]
	write.table(present, file="1_calls/affy/present_mas5.txt", sep="\t", quote=F, col.names=F, row.names=F)
	
	write.table(cbind(names(control.calls),control.calls), file="1_calls/affy/control_calls_mas5.txt", sep="\t", quote=F, col.names=F, row.names=F)
	write.table(cbind(names(het.calls),het.calls), file="1_calls/affy/treatment_calls_mas5.txt", sep="\t", quote=F, col.names=F, row.names=F)
	result<-present
}


##################################################################################################
###### present_gcrma: implementation as described in  
###### 'Schuster et al. "Correcting for sequence biases in present/absent calls." Genome Biology'
##################################################################################################

 
present_gcrma<-function(data) {
	
all.ps <- geneNames(data)
data -> data.trans

# Get the locations the probes in each probesets on the array
pms <- unlist(indexProbes(data, "pm"))
mms <-  unlist(indexProbes(data, "mm"))

pms_test <- indexProbes(data, "pm") #list the probes according to their index - give the lines of the probes on the affybatch
#looks like that:
#> head(pms_test)
#$`1616608_a_at`
# [1] 104972 233841 402186 185398 282536 220898 343563 103877  18388 261917
#[11] 475257  43295   4121 199087
#
#$`1622892_s_at`
# [1] 324552 217447 381823 518150 408044 136102 409879 179444 157504 401895
#[11] 165071 264361 419748  81566
#... - the same can be done for mismatches

#  calculate the single nucleotide affinities that are specific to each replicate sample
ai <- compute.affinities.local(data, Array=NULL)

#  correct for optical background. First do optimal step (setting minimal values for PM/MM). Second transform PM/MM values cel file by cel file.
# in the models for 
op <- bg.adjust.optical(data)

#gc-RMA
#Corrects for background noise, as well as non-specific binding 
#Probe affinity modeled as sum of position-dependent base effects 
#Calculated for each probe based on sequence information 
#PM’s are adjusted by subtracting “shrunken”MM, corrected for its affinity


#  to calculate gc-RMA transformed PM and MM probe values
for(i in 1:length(sampleNames(data))) {
	# what factor do I need the ai values multiply with to get the op values?
	parameters <- bg.parameters.ns(intensity(op)[mms,i], intensity(ai)[mms,i], intensity(ai)[c(pms, mms),i])
	mu.pm <- parameters$bg.mu2
	sigma <- parameters$bg.sigma
	intensity(data.trans)[c(pms, mms),i] <- gcrma.bg.transformation(intensity(op)[c(pms, mms),i], mu.pm, sigma, k = .5)
}
all.exprsVal=expresso(data.trans, bg.correct=F, normalize=F, pmcorrect.method="pmonly",summary.method="avgdiff")

save(data.trans, file='1_calls/data.trans.dat')
save(all.exprsVal, file='1_calls/all.exprsVal.dat')
#  calculate present/absent calls
#  mas5calls Performs the Wilcoxon signed rank-based gene expression presence/absence detection algorithm
#  alpha values set the "present" and "marginal" P value cutoffs
#  alpha1 is set to 0.03 and the expected false discover rate at this P value is roughly 5% to 8%
#  alpha2 is set to 0.12 and the expected false discover rate at this P value is roughly 10% to 13%

trans.calls <- mas5calls(data.trans, tau=0.1, alpha1 = 0.037, alpha2 = 0.111)
#create an expressionSet object of detection values - "P", "M" or "A"
#CONTROL.ARRAYS=c("w1h.CEL","w2h.CEL","w3h.CEL")
#TREATMENT.ARRAYS=c("d1h.CEL","d2h.CEL","d3h.CEL")
CONTROL.ARRAYSTrans<-assayData(trans.calls)[["se.exprs"]][,CONTROL.ARRAYS] #extract standard errors of the ctrl arrays
ablatedArraysTrans<-assayData(trans.calls)[["se.exprs"]][,TREATMENT.ARRAYS] # the same for the treatment
control.trans.calls <- apply(CONTROL.ARRAYSTrans, 1, function(x) median(x)) #calculate median for each line
ablated.trans.calls <- apply(ablatedArraysTrans, 1, function(x) median(x)) # median calculations
present.trans <- all.ps[control.trans.calls < .111 | ablated.trans.calls < .111]# combine the probe names of both lists of present calls (present and marginal) together

CONTROL.ARRAYSnames=all.ps[control.trans.calls < .111]
ablatedArraysnames= all.ps[ablated.trans.calls < .111]
save(present.trans, file='1_calls/present.trans.dat') # save the object of all presents gene ids
save(CONTROL.ARRAYSnames, file='3_simplePierreAlchemy/fractions/controlArraysNames.dat') #save the object of the genes with their present-p-values for the control experiments with a p-val<0.111
save(ablatedArraysnames, file='3_simplePierreAlchemy/fractions/ablatedArraysNames.dat')#save the object of the genes with their present-p-values for the treatment experiments with a p-val<0.111
presenttrans2=apply(assayData(trans.calls) [['se.exprs']],2, function(x) sum(x<.111)) # how many genes I have within these parameters

present.exprsVal=all.exprsVal[present.trans,] # saving the subset of the affybatch for only those probes with a present call.
save(present.exprsVal, file="1_calls/present.exprsVal.dat")
write.table(exprs(present.exprsVal), file="1_calls/expressionValues_presentProbes.txt", sep="\t", quote=F)

length(present.trans)
#cat(length(present.trans));cat(paste(' present.transL', '\n', sep = ''))
#  cat(length(presenttrans2));cat(paste(' presenttrans2L', '\n', sep = ''))
cat(paste('number of genes which passes the threshold for each conditions in present2trans are:', '\n', sep = ''));print(presenttrans2)

#extracting the names of present probe IDs for treatment, control and all together
write.table(ablatedArraysnames, file="1_calls/affy/present_treatment_gcrma.txt", sep="\t", quote=F, col.names=F, row.names=F)
write.table(CONTROL.ARRAYSnames, file="1_calls/affy/present__control_gcrma.txt", sep="\t", quote=F, col.names=F, row.names=F)
write.table(present.trans, file="1_calls/affy/present_gcrma.txt", sep="\t", quote=F, col.names=F, row.names=F)

cat("writing table to 1_calls/affy/present_gcrma.txt"); cat("\n")
write.table(cbind(names(control.trans.calls),control.trans.calls), file="1_calls/affy/control_calls_gcrma.txt", sep="\t", quote=F, col.names=F, row.names=F)
cat("writing table to 1_calls/affy/control_calls_gcrma.txt"); cat("\n")
write.table(cbind(names(ablated.trans.calls),ablated.trans.calls), file="1_calls/affy/treatment_calls_gcrma.txt", sep="\t", quote=F, col.names=F, row.names=F)
cat("writing table to 1_calls/affy/treatment_calls_gcrma.txt}"); cat("\n")
}



library(affy)#environment for data analysis + exploration of Aﬀymetrix oligonucleotide array probe level data 
library(gcrma)

source("0_init.R")
baseDir <-getwd()
all.ps <- geneNames(data) #list of the probe IDs
#present_mas5 <- present_mas5(data)
present_gcrma <- present_gcrma(data)



################################ Highlight overall density of selected files (ex outliers) in graph ############# 
load('1_calls/data.trans.dat')
load('1_calls/present.trans.dat')
cat("loaded files from 1_calls/data.trans.dat and 1_calls/present.trans.dat \n")
out=paste(BASEDIR, '/1_calls/Densities_celHighlight.pdf',sep='')
cat("save pdf file under", out, "\n")
pdf(out)
par(mfrow=c(1,1))

plot=plotDensity(log2(pm(data.trans, present.trans)),  lty=1, cex=14)
dev.off()
print("finished running the calls.R script")



