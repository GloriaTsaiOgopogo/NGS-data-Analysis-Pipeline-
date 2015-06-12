present_gcrma<-function(data) {
	
	all.ps <- geneNames(data)
	data -> data.trans

	# Get the locations the probes in each probesets on the array
	pms <- unlist(indexProbes(data, "pm"))
	mms <-  unlist(indexProbes(data, "mm"))

	pms_test <- indexProbes(data, "pm") 
	#  calculate the single nucleotide affinities that are specific to each replicate sample
	ai <- compute.affinities.local(data, Array=NULL)

	#  correct for optical background. First do optimal step (setting minimal values for PM/MM). Second transform PM/MM values cel file by cel file.
	# in the models for 
	op <- bg.adjust.optical(data)

	#gc-RMA
	#Corrects for background noise, as well as non-specific binding 
	#Probe affinity modeled as sum of position-dependent base effects 
	#Calculated for each probe based on sequence information 
	#PM’s are adjusted by subtracting “shrunken”MM, 
	#corrected for its affinity

	#  to calculate gc-RMA transformed PM and MM probe values
	for(i in 1:length(sampleNames(data))) {
		# what factor do I need the ai values multiply with to get the op values?
		parameters <- bg.parameters.ns(intensity(op)[mms,i], intensity(ai)[mms,i], intensity(ai)[c(pms, mms),i])
		mu.pm <- parameters$bg.mu2
		sigma <- parameters$bg.sigma
		intensity(data.trans)[c(pms, mms),i] <- gcrma.bg.transformation(intensity(op)[c(pms, mms),i], mu.pm, sigma, k = .5)
	}

	#  calculate present/absent calls
	#  alpha values set the "present" and "marginal" P value cutoffs
	#  alpha1 is set to 0.03 and the expected false discover rate at this P value is roughly 5% to 8%
	#  alpha1 is set to 0.12 and the expected false discover rate at this P value is roughly 10% to 13%

	trans.calls <- mas5calls(data.trans, tau=0.1, alpha1 = 0.037, alpha2 = 0.111)

	HIGH.CTR.ARRAYS.Trans<-assayData(trans.calls)[["se.exprs"]][,HIGH.CTR.ARRAYS]
	HIGH.LIT.ARRAYS.Trans<-assayData(trans.calls)[["se.exprs"]][,HIGH.LIT.ARRAYS]

	LOW.CTR.ARRAYS.Trans<-assayData(trans.calls)[["se.exprs"]][,LOW.CTR.ARRAYS]
	LOW.LIT.ARRAYS.Trans<-assayData(trans.calls)[["se.exprs"]][,LOW.LIT.ARRAYS]

	TOTAL.CTR.ARRAYS.Trans<-assayData(trans.calls)[["se.exprs"]][,TOTAL.CTR.ARRAYS]
	TOTAL.LIT.ARRAYS.Trans<-assayData(trans.calls)[["se.exprs"]][,TOTAL.LIT.ARRAYS]

	HIGH.CTR.ARRAYS.calls <- apply(HIGH.CTR.ARRAYS.Trans, 1, function(x) median(x))
	HIGH.LIT.ARRAYS.calls <- apply(HIGH.LIT.ARRAYS.Trans, 1, function(x) median(x))
	LOW.CTR.ARRAYS.calls <- apply(LOW.CTR.ARRAYS.Trans, 1, function(x) median(x))
	LOW.LIT.ARRAYS.calls <- apply(LOW.LIT.ARRAYS.Trans, 1, function(x) median(x))
	TOTAL.CTR.ARRAYS.calls <- apply(TOTAL.CTR.ARRAYS.Trans, 1, function(x) median(x))
	TOTAL.LIT.ARRAYS.calls <- apply(TOTAL.LIT.ARRAYS.Trans, 1, function(x) median(x))

	present.trans <- all.ps[HIGH.CTR.ARRAYS.calls < .111 | HIGH.LIT.ARRAYS.calls < .111 |
	LOW.CTR.ARRAYS.calls < .111 | LOW.LIT.ARRAYS.calls < .111 |
	TOTAL.CTR.ARRAYS.calls < .111 | TOTAL.LIT.ARRAYS.calls < .111]

	####################### SAVE BACKGROUND CORRECTED VALUES ######################## 
	#data.trans.sum=rma(data.trans, normalize=F, background=F)
	#write.table(exprs(data.trans.sum), file="1_bgcorrection/affy/UNION_ctrF-ctrM-highF-highM-lowF-lowM_gcrma.txt", sep="\t", quote=F, col.names=T, row.names=T)

	#write.table(exprs(data.trans.sum)[,CTR.ARRAYS], file="1_bgcorrection/affy/COMBINED_ctr_female_male_gcrma_COND8.txt", sep="\t", quote=F, col.names=T, row.names=T)
	#write.table(exprs(data.trans.sum)[,CTR.FEM.ARRAYS], file="1_bgcorrection/affy/ctr_female_gcrma_COND4.txt", sep="\t", quote=F, col.names=T, row.names=T)
	#write.table(exprs(data.trans.sum)[,CTR.MALE.ARRAYS], file="1_bgcorrection/affy/ctr_male_gcrma_COND4.txt", sep="\t", quote=F, col.names=T, row.names=T)
	
	#write.table(exprs(data.trans.sum)[,HIGH.ARRAYS], file="1_bgcorrection/affy/COMBINED_high_female_male_gcrma.txt_COND8", sep="\t", quote=F, col.names=T, row.names=T)
	#write.table(exprs(data.trans.sum)[,HIGH.FEM.ARRAYS], file="1_bgcorrection/affy/high_female_gcrma_COND4.txt", sep="\t", quote=F, col.names=T, row.names=T)
	#write.table(exprs(data.trans.sum)[,HIGH.MALE.ARRAYS], file="1_bgcorrection/affy/high_male_gcrma_COND4.txt", sep="\t", quote=F, col.names=T, row.names=T)
	
	#write.table(exprs(data.trans.sum)[,LOW.ARRAYS], file="1_bgcorrection/affy/COMBINED_low_female_male_gcrma_COND8.txt", sep="\t", quote=F, col.names=T, row.names=T)
	#write.table(exprs(data.trans.sum)[,LOW.FEM.ARRAYS], file="1_bgcorrection/affy/low_female_gcrma_COND4.txt", sep="\t", quote=F, col.names=T, row.names=T)
	#write.table(exprs(data.trans.sum)[,LOW.MALE.ARRAYS], file="1_bgcorrection/affy/low_male_gcrma_COND4.txt", sep="\t", quote=F, col.names=T, row.names=T)

	######################## SAVE RESULTS FROM P/A ######################## 

	write.table(present.trans, file="1_calls/affy/UNION_totalCtr-totalLith-highCtr-highLith-lowCtr-lowLith_gcrma.txt", sep="\t", quote=F, col.names=F, row.names=F)
	
	write.table(all.ps[HIGH.CTR.ARRAYS.calls  < .111], file="1_calls/affy/ctr_high_gcrma_COND3.txt", sep="\t", quote=F, col.names=F, row.names=F)
	write.table(all.ps[HIGH.LIT.ARRAYS.calls < .111], file="1_calls/affy/lith_high_gcrma_COND3.txt", sep="\t", quote=F, col.names=F, row.names=F)

	write.table(all.ps[LOW.CTR.ARRAYS.calls < .111], file="1_calls/affy/ctr_low_gcrma_COND3.txt", sep="\t", quote=F, col.names=F, row.names=F)
	write.table(all.ps[LOW.LIT.ARRAYS.calls < .111], file="1_calls/affy/lith_low_gcrma_COND3.txt", sep="\t", quote=F, col.names=F, row.names=F)

	write.table(all.ps[TOTAL.CTR.ARRAYS.calls < .111], file="1_calls/affy/ctr_total_gcrma_COND3.txt", sep="\t", quote=F, col.names=F, row.names=F)
	write.table(all.ps[TOTAL.LIT.ARRAYS.calls < .111], file="1_calls/affy/lith_total_gcrma_COND3.txt", sep="\t", quote=F, col.names=F, row.names=F)
}


library(affy)#environment for data analysis + exploration of Aﬀymetrix oligonucleotide array probe level data 
library(gcrma)


path_input_data<-"/Volumes/nfs_research2_remote/projects/collaboration/fly/alzheimersDisease/Luke_translation/AGE_Tain_151211 BLUE/"
filename_phenodata<-"wildtype_lithium.phenodata"
setwd(path_input_data)
###################################################
### 2:  Preprocessing (Loading ets)
###################################################

phenodata <- read.AnnotatedDataFrame(paste(path_input_data,filename_phenodata, sep=""),sep=" ")
data <- ReadAffy(
	filenames=phenodata$filenames,
	sampleNames=sampleNames(phenodata),
	phenoData=phenodata,
	celfile.path=paste(path_input_data, "0_raw/", sep=""))

HIGH.CTR.ARRAYS<-rownames(pData(data)[pData(data)[,3] %in% "control_highTranslation", ])
HIGH.LIT.ARRAYS<-rownames(pData(data)[pData(data)[,3] %in% "lithium_highTranslation", ])

LOW.CTR.ARRAYS<-rownames(pData(data)[pData(data)[,3] %in% "control_lowTranslation", ])
LOW.LIT.ARRAYS<-rownames(pData(data)[pData(data)[,3] %in% "lithium_lowTranslation", ])

TOTAL.CTR.ARRAYS<-rownames(pData(data)[pData(data)[,3] %in% "control_totalRNA",])
TOTAL.LIT.ARRAYS<-rownames(pData(data)[pData(data)[,3] %in% "lithium_totalRNA",])



all.ps <- geneNames(data)
present_gcrma <- present_gcrma(data)



