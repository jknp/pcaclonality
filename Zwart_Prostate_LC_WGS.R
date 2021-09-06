### Zwart Prostate samples LC-WGS primary and metastasis scratchbook

setwd("D:\\D Projects\\201707-02 Primary PCa vs Lymph node metastases\\OKData")

library(QDNAseq)
bins <- getBinAnnotations(binSize=30)
readCounts <- binReadCounts(bins)

readCountsFiltered<- applyFilters(readCounts, residual=T, blacklist=T, chromosomes=NA)
readCountsFiltered <- estimateCorrection(readCountsFiltered)

noisePlot(readCountsFiltered)
copyNumbers<-correctBins(readCountsFiltered) 

copyNumbersNormalized<-normalizeBins(copyNumbers)
copyNumbersSmooth<-smoothOutlierBins(copyNumbersNormalized)

plot(copyNumbersSmooth)

exportBins(copyNumbersSmooth, file=".\\copyNumbers_Prostate_30Kb.txt")
exportBins(copyNumbersSmooth, file=".\\Prostate_30Kb.igv", format="igv")

##### CGHcall analysis 

setwd("D:\\D Projects\\201707-02 Primary PCa vs Lymph node metastases\\OKData")
Data<-read.table(file="copyNumbers_Prostate_30Kb.txt", sep="\t", header=T)

library(CGHcall)
raw <- make_cghRaw(Data)
prep <- preprocess(raw, maxmiss = 0, nchrom = 22) 
nor <-  normalize(prep,method = "median", smoothOutliers = TRUE)  

seg <-  segmentData(nor, method = "DNAcopy",nperm=2000,undo.splits="sdundo",min.width=5,undo.SD=2, clen=25, relSDlong=5)

segnorm <- postsegnormalize(seg,inter=c(-0.4,0.4))
listcalls <- CGHcall(segnorm,nclass=5,robustsig="yes",cellularity=1,ncpus=12)
calls <- ExpandCGHcall(listcalls,segnorm, divide=5, memeff=FALSE)
save(calls, file="Zwart_Prostate.Rdata")  



###### Heatmaps per patient

setwd("D:\\D Projects\\201707-02 Primary PCa vs Lymph node metastases\\OKData")
load("Zwart_Prostate.Rdata")  

Sample_ann<-read.table(file="CFMPB293_key.txt", sep="\t", header=T, stringsAsFactors = F)
Index_samples<-match(paste0("X",Sample_ann[,1]), colnames(calls))
calls<-calls[,Index_samples]

## Double check to see if sorted correctly
table(colnames(calls)==paste0("X",Sample_ann[,1]))

for (i in 1:length(unique(Sample_ann$Patient))){

	Sample<-i
	IndexSelectedSamples<-which(Sample_ann[,4]==Sample)

	pdf(file=paste0("Patient_",i,".pdf"), width=7, height=4)
	source(".\\Additional_Scripts\\overviewPlot_MIN-MAX.R")
	overviewPlot(calls[,IndexSelectedSamples], scaling=0.6)
	dev.off()
}

####### Zoomed in heatmaps #########

source(".\\Additional_Scripts\\overviewPlot_MIN-MAX.R")
#Select chromosome, start bp, end bp
Sel <- c(10, 89600000, 89800000)

#Zoom
tiff(file=paste0("Chrom_", Sel[1], ".tiff"), width=4*480, height=2.6*480)
overviewPlot(calls[which(chromosomes(calls)== Sel[1] & bpstart(calls)< Sel[3] & bpend(calls) > Sel[2])], scaling=1) +
  title(main = paste0("Chromosome ", Sel[1], " ", (Sel[2]/1000), " kb to ", (Sel[3]/1000), " kb"))
dev.off()

#Complete Chrom M only
tiff(file=paste0("P_Chrom_", Sel[1], ".tiff"), width=2*480, height=1.3*480)
#overviewPlot(calls[which(chromosomes(calls)== Sel[1])], scaling=1) + title(main = paste0("Chromosome ", Sel[1]))
overviewPlot(calls[,Ind][which(chromosomes(calls)== Sel[1] & bpstart(calls)< Sel[3] & bpend(calls) > Sel[2])], scaling=1) +
  title(main = paste0("Chromosome ", Sel[1], " ", (Sel[2]/1000), " kb to ", (Sel[3]/1000), " kb"))
dev.off()

## Generate .png file for each profile
Path<-getwd()
source(".\\Additional_Scripts\\Plot_chromosome_lines.R")
setwd(paste(Path))
system(paste("mkdir Segmented"))
setwd(paste0(Path, "\\Segmented\\"))

for (i in 1:ncol(segmented(calls))){
	png(file=paste(Sample_ann[i,3], "_Segmented_", colnames(segmented(calls))[i],".png", sep=""), width=4*480, height=2*480)
	plot(copynumber(calls[,i]), pch=".", cex=2, ylim=c(-2,5), xaxt="n")
	points(segmented(calls[,i]), pch=".", cex=3, col="red")
	Chr.names(chromosomes(calls))
	dev.off()
}
setwd(paste(Path))

##### Overview Heatmaps, Ann_V2##### 
Sample_ann_V2 <- read.table(file="CFMPB293_keyV2.txt", sep="\t", header=T, stringsAsFactors = F)
NVcalls <-calls[,which(Sample_ann_V2$Neo =="0")]

overviewPlot(NVcalls, scaling =1)

Sel <- c(23, 66800000, 67000000)
overviewPlot(NVcalls[which(chromosomes(NVcalls)== Sel[1] & bpstart(NVcalls)< Sel[3] & bpend(NVcalls) > Sel[2])], scaling=1) +
  title(main = paste0("Chromosome ", Sel[1], " ", (Sel[2]/1000), " kb to ", (Sel[3]/1000), " kb"))



source(".\\Additional_Scripts\\overviewPlot_MIN-MAX_Xchromosome_correct.R")
tiff(file="OverviewPlot.tiff", width=4*480, height=2.6*480)
overviewPlot(calls, scaling=1)
dev.off()

## Generate igv readable file
library(focalCall)
source(".\\Additional_Scripts\\igvFiles.R")

igvFiles(calls)
system("FrequencyPlot.igv ")


###################### Additional plots ############################# 
#Frequency plots

source(".\\Additional_Scripts\\FrequencyPlot_Xcor.R")

pdf(file="FrequencyPlots_All_samples.pdf", width=14, height=7)
FrequencyPlot(calls, header="All Sample")
dev.off()

pdf(file="FrequencyPlots_P_and_M_samples_Xcor.pdf", width=14, height=14)
par(mfrow=c(2,1))
FrequencyPlot(calls[,which(Sample_ann[,5]=="P")], header="Primary Sample")
FrequencyPlot(calls[,which(Sample_ann[,5]=="M")], header="Metastatic Sample")
dev.off()

# Heatmap of calls data
source("overviewPlot_calls.R")
overviewPlot_calls(calls)


## Number of segments
getSegments<-function(Segs, Chromo){
  # Count segments
  Segment.count<-matrix(data=0, ncol=ncol(Segs), nrow=nrow(Segs))
  for (j in 1:ncol(Segs)){
    k=1 
    test.data<-Segs[,j]
    for (i in 1:length(test.data)){
      ifelse(test.data[i-1]==test.data[i] & Chromo[i-1]==Chromo[i], k<-k, k<-k+1)
      Segment.count[i,j]<-k
    }
    cat("Sample", j, "finished", "\n")
  }
  return(Segment.count)
}

Segment_number_tmp<-getSegments(segmented(calls), chromosomes(calls))
Segment_number<-Segment_number_tmp[nrow(Segment_number_tmp),]

##########Tumor percentages################
tumorp <- read.table(file="Results_CFMPB293.txt", sep="\t", header=T, stringsAsFactors = F)
index_tumorp <- match(Sample_ann$CF_No,  tumorp$CF.sample)
percentages <- subset(tumorp, select=c(1, 8, 20:21))
percentages <- percentages[index_tumorp,]
percentages$Pair_No <- Sample_ann$Pair_No


#########Clonality analysis################

setwd("D:\\D Projects\\201707-02 Primary PCa vs Lymph node metastases\\OKData")
load("Zwart_Prostate.Rdata")  

library(CGHcall)
library(Clonality)


# Order calls according to patient ID
setwd("D:\\D Projects\\201707-02 Primary PCa vs Lymph node metastases\\OKData")
load("Zwart_Prostate.Rdata")  

Sample_ann<-read.table(file="CFMPB293_key.txt", sep="\t", header=T)
Index_samples<-match(paste0("X",Sample_ann[,1]), colnames(calls))
calls<-calls[,Index_samples]

dataCNA<-CNA(copynumber(calls),chrom=chromosomes(calls),maploc=bpstart(calls),sampleid=colnames(calls))

# Decrease resolution (See documentation), max is 15.000 datapoints
#dataAve<- ave.adj.probes(dataCNA,10)
#dataAve<- ave.adj.probes(dataCNA,20)
dataAve<- ave.adj.probes(dataCNA,40)

ptlist<-Sample_ann$Patient
res40<-clonality.analysis(dataAve, ptlist,  pfreq = NULL, refdata = NULL, nmad = 1,  reference = TRUE, allpairs = TRUE)
histogramPlot(LNsort10[,4], results$refLR[,4])

## Selecting Metastases pairs only from Clonality data

#Find index of LNs, initialize for subsetting
for(i in 1:length(Sample_ann$Tumor)){
  Index_LN <- which(Sample_ann$Tumor == "M")
}
res_df10 <- results$LR
ref_df10 <- results$refLR
LNsort10 <- res_df10[FALSE,]

#Subset for Metastases for LR set, then swap GCF_ID for Pair_No in table
for(i in 1:length(Index_LN)){
  LNsort10 <- rbind(LNsort10, subset(res_df10, res_df10[,2] == paste0("X",Sample_ann$GCF_ID[(Index_LN[i])])))
}
for(i in 1:2){
  Check <- match(LNsort10[,i], paste0("X", Sample_ann[,1]))
  LNsort10[,i] <- Sample_ann$Pair_No[Check]
}

#Plotting for sorted or unsorted
Name <- "Ave_40_LNsorted"
png(file= paste0(Name, ".png"), width=2*480, height=1.3*480)

histogramPlot(LNtab40[,4], ref_df40[,4]) + title(main = Name, line = +0.2)
dev.off()

Name <- "Ave_40_unsorted"
png(file= paste0(Name, ".png"), width=2*480, height=1.3*480)
histogramPlot(res_df[,4], ref_df[,4]) + title(main = Name, line = +0.2)
dev.off()


#Adapted Clonality plot
setwd("D:\\D Projects\\201707-02 Primary PCa vs Lymph node metastases\\OKData")
source(".\\Additional_Scripts\\histogramPlot2.R")
LNtab40 <- LNsort40[which(LNsort40[,10] < 0.01),]
LNtab40ns <- LNsort40[which(LNsort40[,10] >= 0.01),]

#Saved all the above in this file
load("Clonality_Results_ave_probes_10_20_40.RData")

histogramPlot2(LNsort40, ref_df40)


########### Correlation plot Heatmap3###############
library(heatmap3)

PatientSelect <- 11 #(c(1,2,3,7,12,14,15,19,21,25,28))
IndexSelectedSamples <- numeric(0)

for(i in 1:length(PatientSelect)){
  for (j in 1:length(Sample_annotation$Patient)){
    Selection<-which(Sample_annotation[,4] == PatientSelect[i])
    }
  IndexSelectedSamples <- c(IndexSelectedSamples, Selection)
  }

corpat <- cor(segmented(calls[,IndexSelectedSamples]))

pdf(file= "heatmap_p11.pdf", width=4*480, height=2.6*480, res=200)
heatmap3(corpat, method = "ward.D2" , breaks = seq(0, 1, length.out = 1025), scale = "none",showColDendro = F, labRow = Sample_ann$Pair_No[IndexSelectedSamples], 
         labCol = Sample_ann$Pair_No[IndexSelectedSamples])



dev.off()



#____________________________Old_Code________________________________________________________________________
######### Scatterplot #############
library(readr)
library(ggplot2)
df <-read.table(file="copyNumbers_Prostate_30Kb.txt", sep="\t", header=T, stringsAsFactors = F)

#sort df, select two pairs and bind
Index_samples<-match(paste0("X",Sample_ann[,1]), colnames(df))
df <-df[,Index_samples]

visual1 = data.frame(df[,30], df[,31])
visual1 <- setNames(visual1, c("p_log", "m_log"))

visual2 = data.frame(df[,77], df[,78])
visual2 <- setNames(visual2, c("p_log", "m_log"))

visuals = rbind(visual1, visual2)
visuals$vis = c(rep("patient10",89360), rep("patient21", 89360))

ggplot(visuals, aes(p_log, m_log, col = vis)) + geom_jitter(shape = ".", alpha = 0.8) + xlim(-2,2) + ylim(-2,2) + xlab("Primary log2") +
  ylab("Metastasis log2") + ggtitle("log2ratio scatter") + geom_abline(intercept = 0, slope = 1, colour = "black") +
  geom_smooth(method = lm) + theme_bw() + scale_color_manual(values = c("black", "#4699dd"))


############GeneBreak##############
library(GeneBreak)

breakpoints <- getBreakpoints(calls)
breakpointsFiltered <- bpFilter( breakpoints, filter = "CNA-ass", threshold = 0.1)
data("ens.gene.ann.hg19")

# Steps below take processing time - 15 mins
breakpointsAnnotated <- addGeneAnnotation( breakpointsFiltered, ens.gene.ann.hg19)
breakpointGenes <- bpGenes(breakpointsAnnotated)
breakpointStatistics <- bpStats( breakpointGenes,level = "gene", method = "Gilbert", fdr.threshold = 0.001)
breakpointStatistics <- bpStats(breakpointStatistics, level = "feature", method = "BH" )

for(i in 1:22){
  tiff(file= paste0("chr", i, "_bpfreqplot.tif"), width=4*480, height=2.6*480, res=200)
  bpPlot(breakpointStatistics, plot.chr = , plot.ylim = 25, fdr.threshold = 0.001, add.jitter = T)
  dev.off()
}
