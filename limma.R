
library(affy)#environment for data analysis + exploration of Aﬀymetrix oligonucleotide array probe level data 
library(gcrma)
library(limma)

###################################################
### 1: Abeta Settings (Mouse data from KCL)
###################################################
path_input_data<-"Lukered"
filename_phenodata<-"wildtype_lithium.phenodata"
adjusted="BH"


###################################################
### 1.1: load data
###################################################
setwd(path_input_data)
phenodata <- read.AnnotatedDataFrame(paste(path_input_data,filename_phenodata, sep=""),sep=" ")
data <- ReadAffy(phenoData=phenodata,celfile.path="0_raw")

present<-read.delim("1_calls/affy/UNION_totalCtr-totalLith-highCtr-highLith-lowCtr-lowLith_gcrma.txt", sep="\t")[,1]

present_ctr_high<-read.delim("1_calls/affy/ctr_high_gcrma_COND3.txt", sep="\t")[,1]
present_lith_high<-read.delim("1_calls/affy/lith_high_gcrma_COND3.txt", sep="\t")[,1]
present_high<-unique(c(as.character(present_ctr_high),as.character(present_lith_high)))

present_ctr_low<-read.delim("1_calls/affy/ctr_low_gcrma_COND3.txt", sep="\t")[,1]
present_lith_low<-read.delim("1_calls/affy/lith_low_gcrma_COND3.txt", sep="\t")[,1]
present_low<-unique(c(as.character(present_ctr_low),as.character(present_lith_low)))

present_ctr_total<-read.delim("1_calls/affy/ctr_total_gcrma_COND3.txt", sep="\t")[,1]
present_lith_total<-read.delim("1_calls/affy/lith_total_gcrma_COND3.txt", sep="\t")[,1]
present_total<-unique(c(as.character(present_ctr_total),as.character(present_lith_total)))

pData(data)$sex <- factor(pData(data)$sex)
pData(data)$condition <- factor(pData(data)$condition)

present_ctr<-unique(c(as.character(present_ctr_high),as.character(present_ctr_low)))
present_lith<-unique(c(as.character(present_lith_high),as.character(present_lith_low)))



#########################################################################################################################################################
####### Data quality: remove low_aggression_line1_female2; ?control_line_male1?
#########################################################################################################################################################


library("affyPLM")
dataPLM <- fitPLM(data)

par(mar=c(14,5,3,3))
boxplot(dataPLM, main="NUSE",ylim = c(0.95, 1.10),
  outline = FALSE, col="lightblue",las=3,
   whisklty=0, staplelty=0)

par(mar=c(14,5,3,3))
 Mbox(dataPLM, main="RLE", ylim = c(-0.4, 0.4),
   outline = FALSE, col="mistyrose", las=3,
   whisklty=0, staplelty=0)

##!## NOT USABLE: control_totalRNA.1        37_15122011_(Drosophila_2).CEL        control_totalRNA.1        control_totalRNA


#########################################################################################################################################################
####### keep polysome fractions only (Translation data)
#########################################################################################################################################################
exclude=grep("totalRNA",pData(data)$label)
tdata <- data[, -exclude]

teset <- rma(tdata)
teset<-teset[c(as.vector(present_low), as.vector(present_high)),]

expression_teset<-unique(exprs(teset))
colnames(expression_teset)<-pData(data)[pData(data)[,1] %in% colnames(expression_teset), 2]

# Calculate ratio high/low for each replicate
ratios<-apply(expression_teset, 1, function(x){
	res<-c()
	for (i in 1:6){
		ratio=x[6+i]-x[i]		
		res=c(res, ratio)
	}
	res
})
ratios_swap<-t(ratios)

#ratios_swap<-2^(ratios_swap)
#2^(log2(800)-log2(1000))

design<-model.matrix(~0+as.factor(c(rep("Ratio.Ctr", 3), rep("Ratio.Lith", 3))))
colnames(design)<-c("Ratio.Ctr", "Ratio.Lith")


# compare female high aggression versus female controls
contrast<-makeContrasts(Ratio.Lith-Ratio.Ctr,levels=design)

fit<-lmFit(ratios_swap, design)
fit2<-contrasts.fit(fit, contrast)
fit2<-eBayes(fit2)


mtop.fitness<-topTable(fit2, coef=1, adjust=adjusted, number=nrow(teset))
mtop.fitness[,c(2:6)]<-signif(mtop.fitness[,c(2:6)], digits=3)


expression_teset<-expression_teset[mtop.fitness[,1],]


controlMeans<-as.data.frame(apply(ratios_swap, 1, function(x){
	signif(mean(x[grep("control",colnames(ratios_swap))]),4)	
}))
colnames(controlMeans)[1]<-"controlRatioMean"
lithiumMeans<-as.data.frame(apply(ratios_swap, 1, function(x){
	signif(mean(x[grep("lithium",colnames(ratios_swap))]),4)	
}))
colnames(lithiumMeans)[1]<-"lithiumRatioMean"

mtop.fitness<-mtop.fitness[,-c(3,6)]

mtop.fitness<-merge(mtop.fitness,controlMeans,by.x="ID", by.y="row.names", sort=F)
mtop.fitness<-merge(mtop.fitness,lithiumMeans,by.x="ID", by.y="row.names", sort=F)
mtop.fitness<-merge(mtop.fitness,ratios_swap,by.x="ID", by.y="row.names", sort=F)

write.table(mtop.fitness, paste(path_input_data,"11_lm/ad_Lith-Ctr_log2.txt", sep=""),col.names=T, row.names=F, quote=F, sep="\t")

# non log
mtop.fitness[,2]<-2^mtop.fitness[,2]
mtop.fitness[,5:12]<-2^mtop.fitness[,5:12]
write.table(mtop.fitness, paste(path_input_data,"11_lm/ad_Lith-Ctr_absolute.txt", sep=""),col.names=T, row.names=F, quote=F, sep="\t")

pdf(paste(path_input_data,"11_lm/ad_Lith-Ctr_cutoff.pdf", sep=""))
	plot(1:nrow(mtop.fitness), log2(mtop.fitness[,3]), ylab="(log2 pvalue)")
dev.off()

###############EXTEND OUTPUT FILES#
library(drosophila2.db)
ensembl2affy<-unlist(as.list(drosophila2FLYBASE))
ensembl2affy<-as.data.frame(ensembl2affy[!is.na(ensembl2affy)]);colnames(ensembl2affy)[1]<-"Ensembl.Gene.ID"
geneSymbols <- as.data.frame(unlist(as.list(drosophila2SYMBOL)));colnames(geneSymbols)[1]<-"symbols"


fly_limma_resultsfile<- paste(path_input_data,"11_lm/ad_Lith-Ctr_absolute.txt", sep="")
fly_limma_results<-read.delim(fly_limma_resultsfile, sep="\t")
fly2description_resultsfile<- paste(path_input_data,"11_lm/flybase2description_mart_export.txt", sep="")
fly2description<-read.delim(fly2description_resultsfile, sep="\t")

fly2description[,2]<- unlist(sapply(fly2description[,2], function(x){strsplit(as.character(x), " \\[")[[1]][1]}))


fly_limma_results.ext<-merge(fly_limma_results, ensembl2affy, by.x="ID", by.y="row.names", all.x=T, sort=F)
fly_limma_results.ext<-merge(fly_limma_results.ext, geneSymbols, by.x="ID", by.y="row.names", all.x=T, sort=F)
fly_limma_results.ext<-merge(fly_limma_results.ext, fly2description, by.x="Ensembl.Gene.ID", by.y="Ensembl.Gene.ID", all.x=T, sort=F)


fly_limma_results.ext<-fly_limma_results.ext[order(fly_limma_results.ext[,5]),]
write.table(fly_limma_results.ext, paste(path_input_data,"11_lm/ad_Lith-Ctr_absolute_extended.txt", sep=""),col.names=T, row.names=F, quote=F, sep="\t")


#########################################################################################################################################################
####### keep Transcriptome data only (total RNA expression diff)
#########################################################################################################################################################
include=grep("totalRNA",pData(data)$label)
totalExpressionData <- data[, include]

expressionEset <- rma(totalExpressionData)
expressionEset<-expressionEset[c(as.vector(present_total)),]

expression_total<-unique(exprs(expressionEset))
colnames(expression_total)<-pData(data)[pData(data)[,1] %in% colnames(expression_total), 2]

design<-model.matrix(~0+as.factor(c(rep("Ctr", 3), rep("Lith", 3))))
colnames(design)<-c("Ctr", "Lith")


# compare female high aggression versus female controls
contrast<-makeContrasts(Lith-Ctr,levels=design)

fit<-lmFit(expressionEset, design)
fit2<-contrasts.fit(fit, contrast)
fit2<-eBayes(fit2)


mtop.fitness<-topTable(fit2, coef=1, adjust=adjusted, number=nrow(teset))
mtop.fitness[,c(5,6)]<-format(mtop.fitness[,c(5,6)], digits=2)
mtop.fitness[,c(2:4)]<-round(mtop.fitness[,c(2:4)], digits=3)



controlMeans<-as.data.frame(apply(expression_total, 1, function(x){
	round(mean(x[grep("control",colnames(expression_total))]),2)	
}))
colnames(controlMeans)[1]<-"controlMean"
lithiumMeans<-as.data.frame(apply(expression_total, 1, function(x){
	round(mean(x[grep("lithium",colnames(expression_total))]),2)	
}))
colnames(lithiumMeans)[1]<-"lithiumMean"

mtop.fitness<-mtop.fitness[,-c(4,7)]

mtop.fitness<-merge(mtop.fitness,controlMeans,by.x="ID", by.y="row.names", sort=F)
mtop.fitness<-merge(mtop.fitness,lithiumMeans,by.x="ID", by.y="row.names", sort=F)


write.table(mtop.fitness, paste(path_input_data,"11_lm/totalRNA_Lith-Ctr.txt", sep=""),col.names=T, row.names=F, quote=F, sep="\t")

#########################################################################################################################################################
######### Transcriptome data
library(drosophila2.db)
ensembl2affy<-unlist(as.list(drosophila2FLYBASE))
ensembl2affy<-as.data.frame(ensembl2affy[!is.na(ensembl2affy)]);colnames(ensembl2affy)[1]<-"Ensembl.Gene.ID"
geneSymbols <- as.data.frame(unlist(as.list(drosophila2SYMBOL)));colnames(geneSymbols)[1]<-"symbols"


fly_limma_resultsfile<- paste(path_input_data,"11_lm/totalRNA_Lith-Ctr.txt", sep="")
fly_limma_results<-read.delim(fly_limma_resultsfile, sep="\t")
fly2description_resultsfile<- paste(path_input_data,"11_lm/flybase2description_mart_export.txt", sep="")
fly2description<-read.delim(fly2description_resultsfile, sep="\t")

fly2description[,2]<- unlist(sapply(fly2description[,2], function(x){strsplit(as.character(x), " \\[")[[1]][1]}))


fly_limma_results.ext<-merge(fly_limma_results, ensembl2affy, by.x="ID", by.y="row.names", all.x=T, sort=F)
fly_limma_results.ext<-merge(fly_limma_results.ext, geneSymbols, by.x="ID", by.y="row.names", all.x=T, sort=F)
fly_limma_results.ext<-merge(fly_limma_results.ext, fly2description, by.x="Ensembl.Gene.ID", by.y="Ensembl.Gene.ID", all.x=T, sort=F)


fly_limma_results.ext<-fly_limma_results.ext[order(fly_limma_results.ext[,5]),]

pdf(paste(path_input_data,"11_lm/cutoff_selection.pdf", sep=""))
	plot(1:nrow(fly_limma_results.ext), log2(fly_limma_results.ext[,6]), ylab="log2(adj.pValue)")
	plot(1:nrow(fly_limma_results.ext), log2(fly_limma_results.ext[,6]), xlim=c(1,1000), ylab="log2(adj.pValue)")
dev.off()
	
write.table(fly_limma_results.ext, paste(path_input_data,"11_lm/totalRNA_Lith-Ctr_extended.txt", sep=""),col.names=T, row.names=F, quote=F, sep="\t")


