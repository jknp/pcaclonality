##Clonality analysis
setwd("D:\\D Projects\\201707-02 Primary PCa vs Lymph node metastases\\OKData")
load("Zwart_Prostate.Rdata")  

library(CGHcall)
library(Clonality)

# Order calls according to patient ID
setwd("D:\\D Projects\\201707-02 Primary PCa vs Lymph node metastases\\OKData")
load("Zwart_Prostate.Rdata")  

Sample_annotation<-read.table(file="CFMPB293_key.txt", sep="\t", header=T)
Index_samples<-match(paste0("X",Sample_annotation[,1]), colnames(calls))
calls<-calls[,Index_samples]

dataCNA<-CNA(copynumber(calls),chrom=chromosomes(calls),maploc=bpstart(calls),sampleid=colnames(calls))

## Decrease resolution (See documentation)
## Max is 15.000 datapoints

dataAve<- ave.adj.probes(dataCNA,10)
#dataAve<- ave.adj.probes(dataCNA,20)

ptlist<-Sample_annotation$Patient

results<-clonality.analysis(dataAve, ptlist,  pfreq = NULL, refdata = NULL, nmad = 1,  reference = TRUE, allpairs = TRUE)
histogramPlot(results1$LR[,4], results1$refLR[,4])

## Selecting Metastases pairs only from Clonality data

#Find index of LNs, initialize for subsetting
for(i in 1:length(Sample_annotation$Tumor)){
  Index_LN <- which(Sample_annotation$Tumor == "M")
}
res_df <- results$LR
ref_df <- results$refLR
test <- res_df[FALSE,]

#Subset for Metastases for LR set, ref is OK
for(i in 1:length(Index_LN)){
  ind <- Index_LN[i]
  test <- rbind(test, subset(res_df, res_df[,2] == paste0("X",Sample_annotation$GCF_ID[ind])))
}

#Plotting for sorted or unsorted
Name <- "Ave_10_LNsorted"
png(file= paste0(Name, ".png"), width=2*480, height=1.3*480)
histogramPlot(test[,4], ref_df[,4]) + title(main = Name, line = +0.2)
dev.off()

Name <- "Ave_10_unsorted"
png(file= paste0(Name, ".png"), width=2*480, height=1.3*480)
histogramPlot(res_df[,4], ref_df[,4]) + title(main = Name, line = +0.2)
dev.off()

