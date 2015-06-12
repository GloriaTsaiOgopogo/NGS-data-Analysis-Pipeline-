source("0_init.R")
library("drosophila2.db")


#################################################
########  Predefined Function(s)         ########
#################################################
export.catmap.files <- function(simple.pierre.table, is.decreasing=TRUE, type = type){

   #  make sure that the stat column is numeric
    simple.pierre.table$stat<-as.numeric(as.vector(simple.pierre.table$stat))

	# First order it according to absolute stat value: if there are two rows for the same gene: take the gene with the higher stat value
	# This is irrespective of up/downregulated genes, and is only done to check for significance; it is reordered again in line 28 according to
	# the actual value rather than the abslutote value.
	ordered_table<-as.vector(simple.pierre.table[order(abs(as.numeric(as.vector(simple.pierre.table$stat))), decreasing = TRUE),])

  	fbnames<-unique(as.vector(ordered_table$fb.genes))
  	fbnames <-fbnames[fbnames!="NA"]

	# then only keep the line with first occurence
	tmp<-sapply(fbnames,function(x){
		as.numeric(ordered_table[ordered_table$fb.genes %in% x,colnames(ordered_table) %in% "stat" ])[1]
	})

	ordered_table
	# reorder increasing/decreasing
	result<-fbnames[order(tmp, decreasing=as.logical(is.decreasing))]

	write.table(result, file=paste("4_genelists/geneids/",suffix,"_listfile_", type, ".tsv", sep=""), quote=F, row.names=F, col.names=F, eol="\t")
  	write.table(result, file=paste("4_genelists/geneids/",suffix,"_listfile_", type, ".txt", sep=""), quote=F, row.names=F, col.names=F, eol="\n")
}


#################################################
########  Get Up/Downregulted genes      ########
#################################################

    suffix<-"unpaired"
#    numberColumns2export<-8

  	simple_pierre_alchemy <- read.table(file=paste("3_simplePierreAlchemy","unpaired_result_present.tbl", sep="/"))

	#################################################
	########  Load Present  GeneIds          ########
	#################################################

	# get mapping fb2affy
	# drosophila2FLYBASE is an R object that contains mappings between manufacturer identiﬁers and
	fb2affy<-unlist(as.list(drosophila2FLYBASE)) # list of 18952 IDs as a chracter vector 
	fb2affy<-as.data.frame(fb2affy[!is.na(fb2affy)])# remove NA's
# left with 13708 IDs with a FB ID.
	colnames(fb2affy)<-"fb.genes"
	
	fb2fb<-unique(cbind(as.character(fb2affy[,1]), as.character(fb2affy[,1]))) # remove multiple IDs and keep only FB gene IDs
	write.table(fb2fb, file="1_calls/geneids/fb2fb.txt", sep="\t", quote=F, col.names=F, row.names=F)
	
	#convert to Flybase geneids
	#which affy IDs do I have in fb2affy that are also present in simple_pierre_alchemy?
	present_geneid <- unique(fb2affy[rownames(fb2affy) %in% rownames(simple_pierre_alchemy), 1])
	write.table(present_geneid, file="1_calls/geneids/present_gcrma_geneids.txt", sep="\t", quote=F, col.names=F, row.names=F)

 	  for (i in seq(as.double(MINIMUM),as.double(MAXIMUM),as.double(STEPS))){
		#ignores whether there is something or not above the cutoff, writes empty files in genelists/affyids
	    up <- rownames(simple_pierre_alchemy[simple_pierre_alchemy$mean.fc > 0 & simple_pierre_alchemy$q < i, ])
	    down <- rownames(simple_pierre_alchemy[simple_pierre_alchemy$mean.fc < 0 & simple_pierre_alchemy$q < i, ])

	    write.table(up, file=paste("4_genelists/affy/",suffix,"_upregulated_",i, ".txt", sep=""), sep="\t", quote=F, col.names=F, row.names=F)
	    write.table(down, file=paste("4_genelists/affy/",suffix,"_downregulated_",i, ".txt", sep=""), sep="\t", quote=F, col.names=F, row.names=F)

		#################################################
		########  Convert Affy2Flybase GeneIds   ########
		#################################################
	    up_geneid <- unique(fb2affy[rownames(fb2affy) %in% up, 1])
	    down_geneid <- unique(fb2affy[rownames(fb2affy) %in% down, 1])

		up_geneid_nonRed <- up_geneid[!up_geneid %in% down_geneid]
	    down_geneid_nonRed <-down_geneid[!down_geneid %in% up_geneid] 
		
		
		if (length(up_geneid)>0){
		    write.table(up_geneid, file=paste("4_genelists/geneids/",suffix,"_upregulated_",i, ".txt", sep=""), sep="\t", quote=F, col.names=F, row.names=F)
		    write.table(up_geneid_nonRed, file=paste("4_genelists/geneids/",suffix,"_upregulated_",i, "_nonRed.txt", sep=""), sep="\t", quote=F, col.names=F, row.names=F)
		}else{
			cat("No genes available for ");cat(paste("4_genelists/geneids/",suffix,"_upregulated_",i, ".txt", sep=""));cat("\n")
		}

		if (length(down_geneid)>0){
		    write.table(down_geneid, file=paste("4_genelists/geneids/",suffix,"_downregulated_",i, ".txt", sep=""), sep="\t", quote=F, col.names=F, row.names=F)
		    write.table(down_geneid_nonRed, file=paste("4_genelists/geneids/",suffix,"_downregulated_",i, "_nonRed.txt", sep=""), sep="\t", quote=F, col.names=F, row.names=F)
		}else{
			cat("No genes available for ");cat(paste("4_genelists/geneids/",suffix,"_downregulated_",i, ".txt", sep=""));cat("\n")
		}


	    #################################################
		########  add Flybase genes      ################
		#################################################

		# add another column and initialize
	  	simple_pierre_alchemy_ext <- merge(simple_pierre_alchemy, fb2affy, by.x="row.names", by.y="row.names", all.x=T)

	  	#################################################
	  	########  add additional columns      ################
	  	#################################################


	  	# geneNames <- unlist(as.list(drosophila2GENENAME)) # add for version 2.7
	    geneSymbols <- as.data.frame(unlist(as.list(drosophila2SYMBOL)));colnames(geneSymbols)[1]<-"symbols"
		geneEntrezGeneAcc <- as.data.frame(unlist(as.list(drosophila2ACCNUM)));colnames(geneEntrezGeneAcc)[1]<-"entrezGeneAcc"
		geneFlybaseCGAcc <- as.data.frame(unlist(as.list(drosophila2FLYBASECG)));colnames(geneFlybaseCGAcc)[1]<-"flybaseCGAcc"
		geneName <- as.data.frame(unlist(as.list(drosophila2GENENAME)));colnames(geneName)[1]<-"geneName"

	  	simple_pierre_alchemy_ext <- merge(simple_pierre_alchemy_ext, geneSymbols, by.x="Row.names", by.y="row.names", all.x=T)
	  	simple_pierre_alchemy_ext <- merge(simple_pierre_alchemy_ext, geneEntrezGeneAcc, by.x="Row.names", by.y="row.names", all.x=T)
	  	simple_pierre_alchemy_ext <- merge(simple_pierre_alchemy_ext, geneFlybaseCGAcc, by.x="Row.names", by.y="row.names", all.x=T)
	  	simple_pierre_alchemy_ext <- merge(simple_pierre_alchemy_ext, geneName, by.x="Row.names", by.y="row.names", all.x=T)


        simple_pierre_alchemy_ext<-simple_pierre_alchemy_ext[order(simple_pierre_alchemy_ext$q),]
	  	write.table(simple_pierre_alchemy_ext, file=paste("3_simplePierreAlchemy/", suffix, "_result_present_extended.tbl", sep=""), sep="\t", quote=F, col.names=T, row.names=T)
# q was changd from 0.1 (single integer) to i (to fit the loop varaibles)
        simple_pierre_alchemy_ext_down<-simple_pierre_alchemy_ext[simple_pierre_alchemy_ext$mean.fc<0 &simple_pierre_alchemy_ext$q< i,]
	  	write.table(simple_pierre_alchemy_ext_down, file=paste("3_simplePierreAlchemy/", suffix, "_result_present_extended_down.tbl_", i, sep=""), sep="\t", quote=F, col.names=T, row.names=T)
# q was changd from 0.1 (single integer) to i (to fit the loop varaibles)
        simple_pierre_alchemy_ext_up<-simple_pierre_alchemy_ext[simple_pierre_alchemy_ext$mean.fc>0&simple_pierre_alchemy_ext$q<i,]
	  	write.table(simple_pierre_alchemy_ext_up, file=paste("3_simplePierreAlchemy/", suffix, "_result_present_extended_up.tbl_", i, sep=""), sep="\t", quote=F, col.names=T, row.names=T)


	    #################################################################
	    ######## Export Decreading/Increasing gene list  ################
	    #################################################################
		export.catmap.files(simple.pierre.table=simple_pierre_alchemy_ext, is.decreasing=as.logical(FALSE), type="increasing")
		export.catmap.files(simple.pierre.table=simple_pierre_alchemy_ext, is.decreasing=as.logical(TRUE), type="decreasing")
 }

		#################################################################
		######## Plots  ################
		#################################################################
	  pdf(paste(BASEDIR, "/3_simplePierreAlchemy/",suffix,"_plot.pdf", sep=""), width=8, height=12)
	  numberPlots<-(0.3-0.2)/STEPS+1
	  par(mfrow=c(numberPlots,1))
	  for (i in seq(MINIMUM,MAXIMUM,STEPS)) {
#	    plot(simple_pierre_alchemy[abs(simple_pierre_alchemy[,foldchange]) > 0,c(mean.sig,foldchange)], main="Fold-change", pch=19, cex=0.6)
plot( simple_pierre_alchemy$mean.fc,simple_pierre_alchemy$mean.sig, main="Fold-change", pch=19, cex=0.6)
#	    simple_pierre_alchemy_fold1<-simple_pierre_alchemy[simple_pierre_alchemy[,foldchange] > 0||simple_pierre_alchemy[,foldchange] < 0,]
simple_pierre_alchemy_fold1<-simple_pierre_alchemy[simple_pierre_alchemy$mean.fc > 0||simple_pierre_alchemy$mean.fc < 0,]
#	    points(simple_pierre_alchemy_fold1[simple_pierre_alchemy_fold1[,q_value] < i, c(mean.sig,foldchange)], col="blue", pch=19, cex=0.6)
points(simple_pierre_alchemy_fold1$q < i, simple_pierre_alchemy$mean.sig, col="blue", pch=19, cex=0.6)
#	    legend("topright", legend=c(paste(i,'q', sep='')), fill=c("blue"))
legend("topright", legend=c(paste(i,'q', sep='')), fill=c("blue"))
	  }
	  dev.off()


