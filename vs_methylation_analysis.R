## Analysis of vestibular schwannoma methylation data 
library(minfi)
library(DMRcate)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(limma)
library(ggplot2)
library(ggthemes)
library(ggfortify)
library(survival)
library(survminer)
library(data.table)
library(DAAG)
library(xlsx)
library(glmnet)
library(reshape)
library(broom)

source("https://bioconductor.org/biocLite.R")
biocLite("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")

setwd('/media/liuj/data5/raleigh/vestibular/analysis')
baseDir <- "/media/liuj/data5/raleigh/vestibular/DNA_methylation/data/IDAT FILES"

# load array data
targets <- read.metharray.sheet(baseDir)
RGSet <- read.metharray.exp(targets=targets)
RGSet@annotation = c(array= "IlluminaHumanMethylationEPIC", annotation = "ilm10b2.hg19")
phenodata <- pData(RGSet)

# get metadata
metaData <- read.xlsx("VS_master_liu.xlsx", 1)
row.names(metaData) <- as.character(metaData$array_sample)
targetsSmp <- targets
row.names(targetsSmp) <- as.character(targetsSmp$barcode)

# normalize array data
manifest <- getManifest(RGSet)
MSet.funnorm <- preprocessFunnorm(RGSet)

# filter data based on common SNPs, X/Y chromosome, cross reactive probes, detection p-value
# remove SNPs
mSetSqFlt = dropLociWithSnps(MSet.funnorm)

#remove X/Y chromsome
annEpic = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
keep = !(featureNames(mSetSqFlt) %in% annEpic$Name[annEPIC$chr %in% c("chrX","chrY")])
mSetSqFlt = mSetSqFlt[keep,]

# remove cross reactive probes [based on published data]
xReactiveProbes = read.csv(file="/home/liuj/genomes/Homo_sapiens/methylation/48639-non-specific-probes-Illumina450k.csv", stringsAsFactors=FALSE)
keep = !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)
mSetSqFlt = mSetSqFlt[keep,]

# remove probes that do not meet mean detection p-value < 0.05
detP = detectionP(RGSet)
keep = colMeans(detP) < 0.05

#sample p-values QC plot
barplot(colMeans(detP),las=2,cex.names=0.8,main="Mean detection p-values")

detP = detP[match(featureNames(mSetSqFlt),rownames(detP)),]
keep = rowSums(detP < 0.05) == ncol(mSetSqFlt)
mSetSqFlt.2 = mSetSqFlt[keep,]

#get beta values
bVals = getBeta(mSetSqFlt.2)
mVals = getM(mSetSqFlt.2)

# PCA and variable probes
pca <- prcomp(t(bVals),center=TRUE,scale.=TRUE)
prop<-summary(pca)
prop<-prop$importance[2,]

# first 2 PC's explain 10% or greater variance
topgenes<-data.frame(pca$rotation[,1:2],maxabs=apply(pca$rotation[,1:2],1,function(x) max(abs(x))))
topgenes<-topgenes[order(-topgenes$maxabs),]
topgenes<-data.frame(topgenes,maxpc=apply(topgenes[,1:2],1,function(x) which.max(abs(x))))
topgenes <- topgenes[topgenes$maxabs > 0,]
table(topgenes[1:5000,"maxpc"])

# get variances for betas
bValsVariance <- data.frame(row.names = row.names(bVals), probes = row.names(bVals), variance = apply(bVals,1,var))
bValsVariance <- bValsVariance[order(-bValsVariance$variance),] # ranked variances

# get heatmap with metadata top 2000 most variable probes
bValTopMatrix <- bVals[row.names(bValsVariance[1:2000,]),]
colnames(bValTopMatrix) <- as.character(targetsSmp[colnames(bValTopMatrix),"sample"])

bValTopMatrixAnnot <- rbind(t(metaData[colnames(bValTopMatrix),c("array_sample","sex","age","size","EOR","preop_growth_rate_percentyr","postop_growth_rate_percentyr","mnp_classifier","mnp_score","prior_rt","prior_surg","clinically_aggressive","LF","LFFP")]),
                            bValTopMatrix)

# output matrix with metadata # edit rows 
write.table(bValTopMatrixAnnot,'methylation/vs_bVal_top2000_annotv3.txt', sep = '\t', quote=FALSE, row.names=TRUE, col.names=NA)

# get annotated all probes
bValsAnnot <- rbind(t(metaData[as.character(targetsSmp[colnames(bVals),"sample"]),c("array_sample","sex","age","size","EOR","preop_growth_rate_percentyr","postop_growth_rate_percentyr","mnp_classifier","mnp_score","prior_rt","prior_surg","clinically_aggressive","LF","LFFP")]),
                    round(bVals,4))
write.table(bValsAnnot,'methylation/vs_bVal_annotv3.txt', sep = '\t', quote=FALSE, row.names=TRUE, col.names=NA)

# get probe gene names
probeOI <- as.character(read.table('methylation/go/c2_probes.txt')$V1)
annEpicSub = annEpic[match(probeOI,annEpic$Name),c(1:4,12:19,24:ncol(annEpic))]

bedOut <- data.frame(annEpicSub[,c(1:2)],as.numeric(annEpicSub[,2])+1)
write.table(bedOut,'methylation/go/c2_probes.bed',row.names=FALSE,col.names=FALSE,sep='\t', quote=FALSE)

#output gene names
write.table(as.character(annEpicSub$GencodeCompV12_NAME),'methylation/go/c2_genes.txt',row.names=FALSE,col.names=FALSE,sep='\t', quote=FALSE)

####################################
# Differential methylation probes ##
# Compare primary failure vs control

metaDmp <- metaData
metaDmp <- data.frame(metaDmp, secondary=as.numeric(metaDmp$prior_rt == 1 | metaDmp$prior_surg ==1))
metaDmp$secondary[is.na(metaDmp$secondary)] = 0
row.names(metaDmp) <- as.character(metaDmp$barcode)

dmpLF <- dmpFinder(bVals[,smpPri], pheno = metaDmp[smpPri,"LF"], type="continuous")

# plot fold changes
plot <- ggplot() +
  geom_point(data = dmpLF[dmpLF$qval > 0.05,], aes(x = -beta,y = -log(qval,10)),shape=16,size=0.5,alpha=1,color='black') + 
  geom_point(data = dmpLF[dmpLF$qval < 0.05,], aes(x = -beta,y = -log(qval,10)),shape=16,size=0.5,alpha=1,color='gray50') + 
  geom_vline(aes(xintercept=0),linetype='dashed',size=0.25) +
  labs(title = paste("Primary tumors \nlocal failure vs local control\n", nrow(dmpLF[dmpLF$qval < 0.05,]),"- probes FDR < 0.05"), x="log odds LC:LF", y="-log10(FDR)") +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_rect(colour='black'))
ggsave("methylation/dmp/plots/primary_lf_lc_flipped.png", plot, width = 2.2, height = 2.5, units='in', device='png')
ggsave("methylation/dmp/plots/primary_lf_lc_flipped.pdf", plot, width = 2.2, height = 2.5, useDingbats = FALSE, limitsize=FALSE)

#labels
plot <- ggplot() +
  geom_point(data = dmpLF[dmpLF$qval > 0.05,], aes(x = -beta,y = -log(qval,10)),shape=16,size=0.5,alpha=1,color='black') + 
  geom_point(data = dmpLF[dmpLF$qval < 0.05,], aes(x = -beta,y = -log(qval,10)),shape=16,size=0.5,alpha=1,color='gray50') +
  geom_point(data = dmpLF[c("cg06957878","cg19723715"),], aes(x = -beta,y = -log(qval,10)), shape=1, size = 1) +
  geom_vline(aes(xintercept=0),linetype='dashed',size=0.25) +
  labs(title = paste("Primary tumors \nlocal failure vs local control\n", nrow(dmpLF[dmpLF$qval < 0.05,]),"- probes FDR < 0.05"), x="log odds LC:LF", y="-log10(FDR)") +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_rect(colour='black'))
ggsave("methylation/dmp/plots/primary_lf_lc_flipped_GPER.png", plot, width = 2.2, height = 2.5, units='in', device='png')


# output custom heatmap - probes of interest with gene level info (name and probe coordinate in the rows )
probesOI <- row.names(dmpLF[dmpLF$qval < 0.05,])
annEpicSub = annEpic[match(probesOI,annEpic$Name),c(1:4,12:19,24:ncol(annEpic))]
annSimp <- data.frame(annEpicSub[,c(1,2,24)])

# annotated volcano plot
annEpicSub = annEpic[match(row.names(dmpLF),annEpic$Name),c(1:4,12:19,24:ncol(annEpic))]
annSimp <- data.frame(annEpicSub[,c(1,2,24)])
dmpLFAnnot <- data.frame(annSimp,signif(dmpLF,5))
write.table(dmpLFAnnot,'methylation/dmp/plots/primary_lf_volcano_annotated.txt', sep = '\t', quote=FALSE, row.names=TRUE, col.names=NA)

# assemble annotated matrix
bValPriSubMatrix <- bValsPri[probesOI,]
colnames(bValPriSubMatrix) <- as.character(targetsSmp[colnames(bValPriSubMatrix),"sample"])
filler = matrix(data = NA,ncol = 3,nrow = 13)
colnames(filler) <- colnames(annSimp)
row.names(filler) <- c("array_sample","sex","age","EOR","prior_rt","prior_surg","size","LF","LFFP","clinically_aggressive","mnp_classifier","mnp_score","methylcluster")

bValSubMatrixAnnot <- data.frame(rbind(filler, annSimp[probesOI,]),
                                rbind(t(metaData[colnames(bValPriSubMatrix),c("array_sample","sex","age","EOR","prior_rt","prior_surg","size","LF","LFFP","clinically_aggressive","mnp_classifier","mnp_score","methylclusters")]),
                                      bValPriSubMatrix))
# output 
write.table(bValSubMatrixAnnot,'methylation/dmp/vs_bValPriLF_dmp_05_annot.txt', sep = '\t', quote=FALSE, row.names=TRUE, col.names=NA)

##################################
## Local control survival analysis
survTbl <- metaData

for (i in 1:nrow(survTbl)){
  x = survTbl[i,"prior_rt_surg"]
  if (!is.na(x)){
    survTbl[i,"prior_rt"] = NA
    survTbl[i,"prior_surg"] = NA
  }
}
# paused for now 

# linear regression for LF vs LC #
probesOI <- c("cg06957878","cg19723715")
#probesOI <- c("cg06957878")
bValSub <- bVals[probesOI,]
colnames(bValSub) <- as.character(targetsSmp[colnames(bValSub),"sample"])
bValSubMean <- data.frame(row.names=colnames(bValSub), smp = colnames(bValSub), meanclust = apply(bValSub,2,mean))

# histogram
plot <- ggplot() +
  geom_histogram(data=bValSubMean, aes(bValSubMean$meanclust), breaks=seq(0, 1, by = .075)) +
  #scale_y_continuous(expand = c(0, 0), limits=c(0,10))+
  theme_few() +
  labs(title = "GPER Promoter methylation", x="beta", y="count") +
  theme(axis.text=element_text(size=6,colour='black'),text=element_text(size=6,colour='black'),axis.text.x=element_text(angle=0,hjust=1,vjust=0.5),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_rect(colour='black'))
ggsave("methylation/dmp/plots/GPER_2probe_mean_hist.pdf", plot, width = 2.5, height = 2.5, useDingbats = FALSE, limitsize=FALSE)

# prepare data for linear fit
classifierTblSub <- data.frame(metaData, meanclust = bValSubMean[row.names(metaData),"meanclust"])
classifierTblSub <- classifierTblSub[colnames(bValPriTopMatrix),c("LF","sex","age","size","EOR","postop_growth_rate_percentyr","meanclust","methylclusters")]

# drop NA's 
classifierTblSubNoNA <- classifierTblSub[!is.na(classifierTblSub$EOR),]

# linear regression
fit <- lm(LF ~ sex + age + size + EOR + postop_growth_rate_percentyr + meanclust, data=classifierTblSub)
fit <- lm(LF ~ sex + age + size + EOR + postop_growth_rate_percentyr, data=classifierTblSub) # model with clinical information only 
summary(fit)
cv.lm(classifierTblSubNoNA, fit, m=3)

# linear regression updated 7/16/2019
fit <- lm(LF ~ size + EOR + meanclust, data=classifierTblSub)
tidy(fit)
write.table(tidy(fit), 'methylation/dmp/mvr_results_size_eor_gper.txt', sep = '\t', quote=FALSE, row.names=FALSE, col.names=TRUE)

fit <- lm(LF ~ size + EOR + methylclusters+meanclust, data=classifierTblSub)
tidy(fit)
write.table(tidy(fit), 'methylation/dmp/mvr_results_size_eor_gper_cluster.txt', sep = '\t', quote=FALSE, row.names=FALSE, col.names=TRUE)

############################################################
# differentially methylated probes in secondary vs primary #
dmpSec <- dmpFinder(bVals, pheno = metaDmp[,"secondary"], type="continuous")

# plot fold changes
plot <- ggplot() +
  geom_point(data = dmpSec[dmpSec$qval > 0.01,], aes(x = beta,y = -log(qval,10)),shape=16,size=0.5,alpha=0.5,color='gray50') + 
  geom_point(data = dmpSec[dmpSec$qval < 0.01,], aes(x = beta,y = -log(qval,10)),shape=16,size=0.5,alpha=0.5,color='green4') + 
  geom_vline(aes(xintercept=0),linetype='dashed',size=0.25) +
  labs(title = paste("Secondary vs Primary tumors \n", nrow(dmpSec[dmpSec$qval < 0.01,]),"- probes FDR < 0.01"), x="log odds Secondary:Primary", y="-log10(FDR)") +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_rect(colour='black'))
ggsave("methylation/dmp/plots/secondary_01.png", plot, width = 2.2, height = 2.5, units='in', device='png')

# output custom heatmap - probes of interest with gene level info (name and probe coordinate in the rows ) - get top 500
#probesOI <- row.names(dmpSec[dmpSec$qval < 0.01,])
probesOI <- row.names(dmpSec[1:500,])
annEpicSub = annEpic[match(probesOI,annEpic$Name),c(1:4,12:19,24:ncol(annEpic))]
annSimp <- data.frame(annEpicSub[,c(1,2,24)])

# assemble annotated matrix
bValSubMatrix <- bVals[probesOI,]
colnames(bValSubMatrix) <- as.character(targetsSmp[colnames(bValSubMatrix),"sample"])
filler = matrix(data = NA,ncol = 3,nrow = 13)
colnames(filler) <- colnames(annSimp)
row.names(filler) <- c("array_sample","sex","age","size","EOR","mnp_classifier","mnp_score","prior_rt","prior_surg","secondary","clinically_aggressive","LF","LFFP")
metaDmpSmp <- metaDmp
row.names(metaDmpSmp) <- as.character(metaDmpSmp$array_sample)

bValSubMatrixAnnot <- data.frame(rbind(filler, annSimp[probesOI,]),
                                 rbind(t(metaDmpSmp[colnames(bValSubMatrix),c("array_sample","sex","age","size","EOR","mnp_classifier","mnp_score","prior_rt","prior_surg","secondary","clinically_aggressive","LF","LFFP")]),
                                       bValSubMatrix))
# output 
write.table(bValSubMatrixAnnot,'methylation/dmp/vs_bValSecondary_dmp_top500_annot.txt', sep = '\t', quote=FALSE, row.names=TRUE, col.names=NA)

############################################################
# differentially methylated probes in RT vs surgery #
metaDmpTx <- metaDmp[metaDmp$secondary == 1,]
metaDmpTx <- metaDmpTx[is.na(metaDmpTx$prior_rt_surg),]
metaDmpTx$prior_rt[is.na(metaDmpTx$prior_rt)] = 0

dmpTx <- dmpFinder(bVals[,row.names(metaDmpTx)], pheno = metaDmpTx[,"prior_rt"], type="continuous")

# plot fold changes
plot <- ggplot() +
  geom_point(data = dmpTx[dmpTx$qval > 0.05,], aes(x = beta,y = -log(qval,10)),shape=16,size=0.5,alpha=0.5,color='gray50') + 
  geom_point(data = dmpTx[dmpTx$qval < 0.05,], aes(x = beta,y = -log(qval,10)),shape=16,size=0.5,alpha=0.5,color='green4') + 
  geom_vline(aes(xintercept=0),linetype='dashed',size=0.25) +
  labs(title = paste("Prior RT vs Prior Surg \n", nrow(dmpSec[dmpTx$qval < 0.05,]),"- probes FDR < 0.05"), x="log odds RT:Surg", y="-log10(FDR)") +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_rect(colour='black'))
ggsave("methylation/dmp/plots/RT_vs_surg_05.png", plot, width = 2.2, height = 2.5, units='in', device='png')

# output custom heatmap - probes of interest with gene level info (name and probe coordinate in the rows )
probesOI <- row.names(dmpTx[dmpTx$qval < 0.05,])
annEpicSub = annEpic[match(probesOI,annEpic$Name),c(1:4,12:19,24:ncol(annEpic))]
annSimp <- data.frame(annEpicSub[,c(1,2,24)])

############################################################
# differentially methylated probes in RT vs no RT #
metaDmpTx <- metaDmp
metaDmpTx$prior_rt[is.na(metaDmpTx$prior_rt)] = 0

dmpTx <- dmpFinder(bVals, pheno = metaDmpTx[,"prior_rt"], type="continuous")

# plot fold changes
plot <- ggplot() +
  geom_point(data = dmpTx[dmpTx$qval > 0.05,], aes(x = beta,y = -log(qval,10)),shape=16,size=0.5,alpha=0.5,color='gray50') + 
  geom_point(data = dmpTx[dmpTx$qval < 0.05,], aes(x = beta,y = -log(qval,10)),shape=16,size=0.5,alpha=0.5,color='green4') + 
  geom_vline(aes(xintercept=0),linetype='dashed',size=0.25) +
  labs(title = paste("Prior RT vs no RT \n", nrow(dmpSec[dmpTx$qval < 0.05,]),"- probes FDR < 0.05"), x="log odds RT:noRT", y="-log10(FDR)") +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_rect(colour='black'))
ggsave("methylation/dmp/plots/RT_vs_noRT_05.png", plot, width = 2.2, height = 2.5, units='in', device='png')

# output custom heatmap - probes of interest with gene level info (name and probe coordinate in the rows )
probesOI <- row.names(dmpTx[dmpTx$qval < 0.001,])
annEpicSub = annEpic[match(probesOI,annEpic$Name),c(1:4,12:19,24:ncol(annEpic))]
annSimp <- data.frame(annEpicSub[,c(1,2,24)])

# assemble annotated matrix
bValSubMatrix <- bVals[probesOI,]
colnames(bValSubMatrix) <- as.character(targetsSmp[colnames(bValSubMatrix),"sample"])
filler = matrix(data = NA,ncol = 3,nrow = 15)
colnames(filler) <- colnames(annSimp)
row.names(filler) <- c("array_sample","sex","age","size","EOR","preop_growth_rate_percentyr","postop_growth_rate_percentyr","mnp_classifier","mnp_score","methylclusters","prior_rt","prior_surg","clinically_aggressive","LF","LFFP")

bValSubMatrixAnnot <- data.frame(rbind(filler, annSimp[probesOI,]),
                                 rbind(t(metaData[colnames(bValSubMatrix),c("array_sample","sex","age","size","EOR","preop_growth_rate_percentyr","postop_growth_rate_percentyr","mnp_classifier","mnp_score","methylclusters","prior_rt","prior_surg","clinically_aggressive","LF","LFFP")]),
                                       bValSubMatrix))
# output 
write.table(bValSubMatrixAnnot,'methylation/dmp/vs_bValRT_dmp_001_annot.txt', sep = '\t', quote=FALSE, row.names=TRUE, col.names=NA)

#######################################################
# differentially methylated probes in RT vs untreated #
metaDmpTx <- metaDmp
metaDmpTx$prior_rt[is.na(metaDmpTx$prior_rt)] = 0
metaDmpTx <- metaDmpTx[is.na(metaDmpTx$prior_surg),]

dmpTx <- dmpFinder(bVals[,row.names(metaDmpTx)], pheno = metaDmpTx[,"prior_rt"], type="continuous")

# plot fold changes
plot <- ggplot() +
  geom_point(data = dmpTx[dmpTx$qval > 0.05,], aes(x = beta,y = -log(qval,10)),shape=16,size=0.5,alpha=0.5,color='gray50') + 
  geom_point(data = dmpTx[dmpTx$qval < 0.05,], aes(x = beta,y = -log(qval,10)),shape=16,size=0.5,alpha=0.5,color='green4') + 
  geom_vline(aes(xintercept=0),linetype='dashed',size=0.25) +
  labs(title = paste("Prior RT only vs primary  \n", nrow(dmpSec[dmpTx$qval < 0.05,]),"- probes FDR < 0.05"), x="log odds RT:noRT", y="-log10(FDR)") +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_rect(colour='black'))
ggsave("methylation/dmp/plots/RT_vs_noTX_05.png", plot, width = 2.2, height = 2.5, units='in', device='png')

plot <- ggplot() +
  geom_point(data = dmpTx[dmpTx$qval > 0.001,], aes(x = beta,y = -log(qval,10)),shape=16,size=0.5,alpha=0.5,color='gray50') + 
  geom_point(data = dmpTx[dmpTx$qval < 0.001,], aes(x = beta,y = -log(qval,10)),shape=16,size=0.5,alpha=0.5,color='green4') + 
  geom_vline(aes(xintercept=0),linetype='dashed',size=0.25) +
  labs(title = paste("Prior RT only vs primary  \n", nrow(dmpSec[dmpTx$qval < 0.001,]),"- probes FDR < 0.001"), x="log odds RT:noRT", y="-log10(FDR)") +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_rect(colour='black'))
ggsave("methylation/dmp/plots/RT_vs_noTX_001.png", plot, width = 2.2, height = 2.5, units='in', device='png')

# output custom heatmap - probes of interest with gene level info (name and probe coordinate in the rows )
probesOI <- row.names(dmpTx[dmpTx$qval < 0.05,])
annEpicSub = annEpic[match(probesOI,annEpic$Name),c(1:4,12:19,24:ncol(annEpic))]
annSimp <- data.frame(annEpicSub[,c(1,2,24)])

############################################################
# differentially methylated probes in surgery vs untreated #
metaDmpTx <- metaDmp
metaDmpTx$prior_surg[is.na(metaDmpTx$prior_surg)] = 0
metaDmpTx <- metaDmpTx[is.na(metaDmpTx$prior_rt),]

dmpTx <- dmpFinder(bVals[,row.names(metaDmpTx)], pheno = metaDmpTx[,"prior_surg"], type="continuous")

# plot fold changes
plot <- ggplot() +
  geom_point(data = dmpTx[dmpTx$qval > 0.05,], aes(x = beta,y = -log(qval,10)),shape=16,size=0.5,alpha=0.5,color='gray50') + 
  geom_point(data = dmpTx[dmpTx$qval < 0.05,], aes(x = beta,y = -log(qval,10)),shape=16,size=0.5,alpha=0.5,color='green4') + 
  geom_vline(aes(xintercept=0),linetype='dashed',size=0.25) +
  labs(title = paste("Prior surgery only vs primary  \n", nrow(dmpSec[dmpTx$qval < 0.05,]),"- probes FDR < 0.05"), x="log odds Surg:noSurg", y="-log10(FDR)") +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_rect(colour='black'))
ggsave("methylation/dmp/plots/surg_vs_noTX_05.png", plot, width = 2.2, height = 2.5, units='in', device='png')

plot <- ggplot() +
  geom_point(data = dmpTx[dmpTx$qval > 0.001,], aes(x = beta,y = -log(qval,10)),shape=16,size=0.5,alpha=0.5,color='gray50') + 
  geom_point(data = dmpTx[dmpTx$qval < 0.001,], aes(x = beta,y = -log(qval,10)),shape=16,size=0.5,alpha=0.5,color='green4') + 
  geom_vline(aes(xintercept=0),linetype='dashed',size=0.25) +
  labs(title = paste("Prior surgery only vs primary  \n", nrow(dmpSec[dmpTx$qval < 0.001,]),"- probes FDR < 0.001"), x="log odds Surg:noSurg", y="-log10(FDR)") +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_rect(colour='black'))
ggsave("methylation/dmp/plots/surg_vs_noTX_001.png", plot, width = 2.2, height = 2.5, units='in', device='png')

# output custom heatmap - probes of interest with gene level info (name and probe coordinate in the rows )
probesOI <- row.names(dmpTx[dmpTx$qval < 0.05,])
annEpicSub = annEpic[match(probesOI,annEpic$Name),c(1:4,12:19,24:ncol(annEpic))]
annSimp <- data.frame(annEpicSub[,c(1,2,24)])

############################################################
# differentially methylated probes in male vs female #
metaDmpTx <- data.frame(metaDmp,stringsAsFactors = FALSE)
metaDmpTx$sex <- as.character(metaDmpTx$sex)
metaDmpTx$sex[metaDmpTx$sex == "M"] = 1
metaDmpTx$sex[metaDmpTx$sex == "F"] = 0

dmpTx <- dmpFinder(bVals[,row.names(metaDmpTx)], pheno = metaDmpTx[,"sex"], type="continuous")

# plot fold changes
plot <- ggplot() +
  geom_point(data = dmpTx[dmpTx$qval >= 0.05,], aes(x = beta,y = -log(qval,10)),shape=16,size=0.5,alpha=0.5,color='gray50') + 
  geom_point(data = dmpTx[dmpTx$qval < 0.05,], aes(x = beta,y = -log(qval,10)),shape=16,size=0.5,alpha=0.5,color='green4') + 
  geom_vline(aes(xintercept=0),linetype='dashed',size=0.25) +
  labs(title = paste("Male vs Female  \n", nrow(dmpSec[dmpTx$qval < 0.05,]),"- probes FDR < 0.05"), x="log odds Male:Female", y="-log10(FDR)") +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_rect(colour='black'))
ggsave("methylation/dmp/plots/M_vs_F_05.png", plot, width = 2.2, height = 2.5, units='in', device='png')

# output custom heatmap - probes of interest with gene level info (name and probe coordinate in the rows )
probesOI <- row.names(dmpTx[dmpTx$qval < 0.05,])
annEpicSub = annEpic[match(probesOI,annEpic$Name),c(1:4,12:19,24:ncol(annEpic))]
annSimp <- data.frame(annEpicSub[,c(1,2,24)])

# assemble annotated matrix
bValSubMatrix <- bVals[probesOI,]
colnames(bValSubMatrix) <- as.character(targetsSmp[colnames(bValSubMatrix),"sample"])
filler = matrix(data = NA,ncol = 3,nrow = 14)
colnames(filler) <- colnames(annSimp)
row.names(filler) <- c("array_sample","sex","age","size","EOR","preop_growth_rate_percentyr","postop_growth_rate_percentyr","mnp_classifier","mnp_score","prior_rt","prior_surg","clinically_aggressive","LF","LFFP")

bValSubMatrixAnnot <- data.frame(rbind(filler, annSimp[probesOI,]),
                                 rbind(t(metaData[colnames(bValSubMatrix),c("array_sample","sex","age","size","EOR","preop_growth_rate_percentyr","postop_growth_rate_percentyr","mnp_classifier","mnp_score","prior_rt","prior_surg","clinically_aggressive","LF","LFFP")]),
                                       bValSubMatrix))
# output 
write.table(bValSubMatrixAnnot,'methylation/dmp/vs_bValSex_dmp_05_annot.txt', sep = '\t', quote=FALSE, row.names=TRUE, col.names=NA)

############################################################
# differentially methylated probes in old vs young #
metaDmpTx <- data.frame(metaDmp,stringsAsFactors = FALSE)
metaDmpTx$age <- as.numeric(metaDmpTx$age)
metaDmpTx$age[metaDmpTx$age < median(metaDmpTx$age)] = 0
metaDmpTx$age[metaDmpTx$age >= median(metaDmpTx$age)] = 1

dmpTx <- dmpFinder(bVals[,row.names(metaDmpTx)], pheno = metaDmpTx[,"age"], type="continuous")

# plot fold changes
plot <- ggplot() +
  geom_point(data = dmpTx[dmpTx$qval >= 0.05,], aes(x = beta,y = -log(qval,10)),shape=16,size=0.5,alpha=0.5,color='gray50') + 
  geom_point(data = dmpTx[dmpTx$qval < 0.05,], aes(x = beta,y = -log(qval,10)),shape=16,size=0.5,alpha=0.5,color='green4') + 
  geom_vline(aes(xintercept=0),linetype='dashed',size=0.25) +
  labs(title = paste("Old vs Young (54.8 yr)  \n", nrow(dmpSec[dmpTx$qval < 0.05,]),"- probes FDR < 0.05"), x="log odds old:young", y="-log10(FDR)") +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_rect(colour='black'))
ggsave("methylation/dmp/plots/old_vs_young_05.png", plot, width = 2.2, height = 2.5, units='in', device='png')

# output custom heatmap - probes of interest with gene level info (name and probe coordinate in the rows )
probesOI <- row.names(dmpTx[dmpTx$qval < 0.05,])
annEpicSub = annEpic[match(probesOI,annEpic$Name),c(1:4,12:19,24:ncol(annEpic))]
annSimp <- data.frame(annEpicSub[,c(1,2,24)])

# assemble annotated matrix
bValSubMatrix <- bVals[probesOI,]
colnames(bValSubMatrix) <- as.character(targetsSmp[colnames(bValSubMatrix),"sample"])
filler = matrix(data = NA,ncol = 3,nrow = 14)
colnames(filler) <- colnames(annSimp)
row.names(filler) <- c("array_sample","sex","age","size","EOR","preop_growth_rate_percentyr","postop_growth_rate_percentyr","mnp_classifier","mnp_score","prior_rt","prior_surg","clinically_aggressive","LF","LFFP")

bValSubMatrixAnnot <- data.frame(rbind(filler, annSimp[probesOI,]),
                                 rbind(t(metaData[colnames(bValSubMatrix),c("array_sample","sex","age","size","EOR","preop_growth_rate_percentyr","postop_growth_rate_percentyr","mnp_classifier","mnp_score","prior_rt","prior_surg","clinically_aggressive","LF","LFFP")]),
                                       bValSubMatrix))
# output 
write.table(bValSubMatrixAnnot,'methylation/dmp/vs_bValAge_dmp_05_annot.txt', sep = '\t', quote=FALSE, row.names=TRUE, col.names=NA)

############################################################
# differentially methylated probes in large vs small #
metaDmpTx <- data.frame(metaDmp,stringsAsFactors = FALSE)
metaDmpTx$size <- as.numeric(metaDmpTx$size)
metaDmpTx$size[metaDmpTx$size < median(metaDmpTx$size,na.rm=TRUE)] = 0
metaDmpTx$size[metaDmpTx$size >= median(metaDmpTx$size,na.rm=TRUE)] = 1
metaDmpTx <- metaDmpTx[!is.na(metaDmpTx$size),]

dmpTx <- dmpFinder(bVals[,row.names(metaDmpTx)], pheno = metaDmpTx[,"size"], type="continuous")

# plot fold changes
plot <- ggplot() +
  geom_point(data = dmpTx[dmpTx$qval >= 0.05,], aes(x = beta,y = -log(qval,10)),shape=16,size=0.5,alpha=0.5,color='gray50') + 
  geom_point(data = dmpTx[dmpTx$qval < 0.05,], aes(x = beta,y = -log(qval,10)),shape=16,size=0.5,alpha=0.5,color='green4') + 
  geom_vline(aes(xintercept=0),linetype='dashed',size=0.25) +
  labs(title = paste("Large vs Small (8.8 cm3)  \n", nrow(dmpSec[dmpTx$qval < 0.05,]),"- probes FDR < 0.05"), x="log odds large:small", y="-log10(FDR)") +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_rect(colour='black'))
ggsave("methylation/dmp/plots/large_vs_small_05.png", plot, width = 2.2, height = 2.5, units='in', device='png')

# output custom heatmap - probes of interest with gene level info (name and probe coordinate in the rows )
probesOI <- row.names(dmpTx[dmpTx$qval < 0.05,])
annEpicSub = annEpic[match(probesOI,annEpic$Name),c(1:4,12:19,24:ncol(annEpic))]
annSimp <- data.frame(annEpicSub[,c(1,2,24)])

# assemble annotated matrix
bValSubMatrix <- bVals[probesOI,]
colnames(bValSubMatrix) <- as.character(targetsSmp[colnames(bValSubMatrix),"sample"])
filler = matrix(data = NA,ncol = 3,nrow = 14)
colnames(filler) <- colnames(annSimp)
row.names(filler) <- c("array_sample","sex","age","size","EOR","preop_growth_rate_percentyr","postop_growth_rate_percentyr","mnp_classifier","mnp_score","prior_rt","prior_surg","clinically_aggressive","LF","LFFP")

bValSubMatrixAnnot <- data.frame(rbind(filler, annSimp[probesOI,]),
                                 rbind(t(metaData[colnames(bValSubMatrix),c("array_sample","sex","age","size","EOR","preop_growth_rate_percentyr","postop_growth_rate_percentyr","mnp_classifier","mnp_score","prior_rt","prior_surg","clinically_aggressive","LF","LFFP")]),
                                       bValSubMatrix))
# output 
write.table(bValSubMatrixAnnot,'methylation/dmp/vs_bValSize_dmp_05_annot.txt', sep = '\t', quote=FALSE, row.names=TRUE, col.names=NA)

####################################################################
# differentially methylated probes in fast vs slow growing post op #
metaDmpTx <- data.frame(metaData,stringsAsFactors = FALSE)
metaDmpTx$postop_growth_rate_percentyr <- as.numeric(metaDmpTx$postop_growth_rate_percentyr)
metaDmpTx$postop_growth_rate_percentyr[metaDmpTx$postop_growth_rate_percentyr <= 0] = 0
metaDmpTx$postop_growth_rate_percentyr[metaDmpTx$postop_growth_rate_percentyr > 0] = 1
metaDmpTx <- metaDmpTx[!is.na(metaDmpTx$postop_growth_rate_percentyr),]

dmpTx <- dmpFinder(bVals[,as.character(metaData[row.names(metaDmpTx),"barcode"])], pheno = metaDmpTx[,"postop_growth_rate_percentyr"], type="continuous")

# plot fold changes
plot <- ggplot() +
  geom_point(data = dmpTx[dmpTx$qval >= 0.05,], aes(x = beta,y = -log(qval,10)),shape=16,size=0.5,alpha=0.5,color='gray50') + 
  geom_point(data = dmpTx[dmpTx$qval < 0.05,], aes(x = beta,y = -log(qval,10)),shape=16,size=0.5,alpha=0.5,color='green4') + 
  geom_vline(aes(xintercept=0),linetype='dashed',size=0.25) +
  labs(title = paste("post op enlarging vs stable/shrink  \n", nrow(dmpSec[dmpTx$qval < 0.05,]),"- probes FDR < 0.05"), x="log odds grow:shrink", y="-log10(FDR)") +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_rect(colour='black'))
ggsave("methylation/dmp/plots/grow_vs_shrink_05.png", plot, width = 2.2, height = 2.5, units='in', device='png')

# output custom heatmap - probes of interest with gene level info (name and probe coordinate in the rows )
probesOI <- row.names(dmpTx[dmpTx$qval < 0.05,])
annEpicSub = annEpic[match(probesOI,annEpic$Name),c(1:4,12:19,24:ncol(annEpic))]
annSimp <- data.frame(annEpicSub[,c(1,2,24)])

# assemble annotated matrix
bValSubMatrix <- bVals[probesOI,]
colnames(bValSubMatrix) <- as.character(targetsSmp[colnames(bValSubMatrix),"sample"])
filler = matrix(data = NA,ncol = 3,nrow = 14)
colnames(filler) <- colnames(annSimp)
row.names(filler) <- c("array_sample","sex","age","size","EOR","preop_growth_rate_percentyr","postop_growth_rate_percentyr","mnp_classifier","mnp_score","prior_rt","prior_surg","clinically_aggressive","LF","LFFP")

bValSubMatrixAnnot <- data.frame(rbind(filler, annSimp[probesOI,]),
                                 rbind(t(metaData[colnames(bValSubMatrix),c("array_sample","sex","age","size","EOR","preop_growth_rate_percentyr","postop_growth_rate_percentyr","mnp_classifier","mnp_score","prior_rt","prior_surg","clinically_aggressive","LF","LFFP")]),
                                       bValSubMatrix))
# output 
write.table(bValSubMatrixAnnot,'methylation/dmp/vs_bValGrowthShrink_dmp_05_annot.txt', sep = '\t', quote=FALSE, row.names=TRUE, col.names=NA)


############################################################################
# output top 2000 gene heatmap with probe annotations for data exploration #
probesOI <- row.names(bValsVariance[1:2000,])
annEpicSub = annEpic[match(probesOI,annEpic$Name),c(1:4,12:19,24:ncol(annEpic))]
annSimp <- data.frame(annEpicSub[,c(1,2,24)])

bValTopMatrix <- bVals[probesOI,]
colnames(bValTopMatrix) <- as.character(targetsSmp[colnames(bValTopMatrix),"sample"])
bValTopMatrixAnnotProbes <- rbind(t(metaData[colnames(bValTopMatrix),c("array_sample","sex","age","size","EOR","preop_growth_rate_percentyr","postop_growth_rate_percentyr","mnp_classifier","mnp_score","prior_rt","prior_surg","clinically_aggressive","LF","LFFP")]),
                            bValTopMatrix)

# assemble annotated matrix with probe data
filler = matrix(data = NA,ncol = 3,nrow = 14)
colnames(filler) <- colnames(annSimp)
row.names(filler) <- c("array_sample","sex","age","size","EOR","preop_growth_rate_percentyr","postop_growth_rate_percentyr","mnp_classifier","mnp_score","prior_rt","prior_surg","clinically_aggressive","LF","LFFP")

bValTopMatrixAnnotProbes <- data.frame(rbind(filler, annSimp[probesOI,]), bValTopMatrixAnnotProbes)

# outpit matrix with metadata # edit rows 
write.table(bValTopMatrixAnnotProbes,'methylation/vs_bVal_top2000_annot_probeannot_v3.txt', sep = '\t', quote=FALSE, row.names=TRUE, col.names=NA)

#####################################
## Analysis of global methylation ###
smpOrder1 <- c('S1_2', 'S1_1', 'S23', 'S38', 'S51', 'S41', 'S31', 'S15', 'S44', 'S49', 'S46', 'S56', 'S16', 'S13', 'S19', 'S28', 'S55', 'S37', 'S50', 'S47', 'S26', 'S22', 'S59_2', 'S21', 'S53', 'S17', 'S66', 'S61_1', 'S61_2', 'S9_1', 'S6_1', 'S52', 'S43', 'S45', 'S60_1', 'S40', 'S3_1', 'S3_2', 'S36', 'S64', 'S39', 'S18', 'S63', 'S61_3', 'S27', 'S30', 'S24', 'S48', 'S54', 'S42', 'S65', 'S58_2', 'S60_2', 'S9_2', 'S29', 'S57', 'S6_2', 'S20', 'S14', 'S67', 'S62', 'S12', 'S11', 'S10', 'S68', 'S69')
bValsPlot <- bVals
colnames(bValsPlot) <- as.character(targetsSmp[colnames(bValsPlot),"sample"])
bValsPlot <- bValsPlot[,smpOrder1]
bValsPlot <- melt(bValsPlot)
colnames(bValsPlot) <- c("probe","smp","b")
bValsPlot$smp <- factor(bValsPlot$smp, levels = smpOrder1)

plot<-ggplot(bValsPlot)+
  geom_boxplot(aes(x=smp,y=b),width=0.9, fill='goldenrod1', outlier.size=0.5, size=0.25, fatten = 1)+
  theme_few() + 
  theme(axis.text=element_text(size=6),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),text=element_text(size=6),axis.ticks=element_line(colour='black',size=0.25),panel.border=element_rect(colour='black')) +
  guides(fill=FALSE)+
  labs(y='beta',x='',title="")
ggsave('methylation/hclust/bval_boxplots.pdf',plot,width=8,height=2.2,useDingbats=FALSE)

#################################################
## Matched pairs analysis of samples ############
smpOrder2 <- c('S1_1', 'S1_2', 'S3_1', 'S3_2', 'S6_1', 'S6_2', 'S9_1', 'S9_2', 'S60_1', 'S60_2', 'S61_1', 'S61_2', 'S61_3')

bValsPaired <- bVals
colnames(bValsPaired) <- as.character(targetsSmp[colnames(bValsPaired),"sample"])
bValsPaired <- bValsPaired[,smpOrder2]

bValsPairedVar <- data.frame(row.names = row.names(bValsPaired), probes = row.names(bValsPaired), variance = apply(bValsPaired,1,var))
bValsPairedVar <- bValsPairedVar[order(-bValsPairedVar$variance),]

bValsPairedTop <- bValsPaired[row.names(bValsVariance[1:2000,]),]
heatmap(cor(bValsPairedTop))
heatmap(cor(bValsPaired))

# output 
probesOI <- row.names(bValsPairedTop)
annEpicSub = annEpic[match(probesOI,annEpic$Name),c(1:4,12:19,24:ncol(annEpic))]
annSimp <- data.frame(annEpicSub[,c(1,2,24)])

bValsPairedTopAnnotProbes <- rbind(t(metaData[colnames(bValsPairedTop),c("array_sample","sex","age","size","EOR","preop_growth_rate_percentyr","postop_growth_rate_percentyr","mnp_classifier","mnp_score","methylclusters","prior_rt","prior_surg","clinically_aggressive","LF","LFFP")]),
                                   bValsPairedTop)

# assemble annotated matrix with probe data
filler = matrix(data = NA,ncol = 3,nrow = 15)
colnames(filler) <- colnames(annSimp)
row.names(filler) <- c("array_sample","sex","age","size","EOR","preop_growth_rate_percentyr","postop_growth_rate_percentyr","mnp_classifier","mnp_score","methylclusters","prior_rt","prior_surg","clinically_aggressive","LF","LFFP")


bValsPairedTopAnnotProbes <- data.frame(rbind(filler, annSimp[probesOI,]), bValsPairedTopAnnotProbes)

# outpit matrix with metadata # edit rows 
write.table(bValsPairedTopAnnotProbes,'methylation/vs_bVal_top2000_pairs_probeannot.txt', sep = '\t', quote=FALSE, row.names=TRUE, col.names=NA)

###############################################
# Pearson distributions of all relevant pairs #
probesOI <- row.names(bValsVariance[1:2000,])
bValTopMatrix <- bVals[probesOI,]
colnames(bValTopMatrix) <- as.character(targetsSmp[colnames(bValTopMatrix),"sample"])

# get primary samples 
priCors <- cor(bValTopMatrix[,colnames(bValPriTopMatrix)])
priCorsTri <- lower.tri(priCors, diag=TRUE)
priCors[priCorsTri] <- NA
priCors <- as.numeric(priCors)
priCors <- priCors[!is.na(priCors)]

# get all prior RT 
rtCors <- cor(bValTopMatrix[,!is.na(as.character(metaData[metaData$prior_rt == 1, "array_sample"]))])
rtCorsTri <- lower.tri(rtCors, diag=TRUE)
rtCors[rtCorsTri] <- NA
rtCors <- as.numeric(rtCors)
rtCors <- rtCors[!is.na(rtCors)]

# get all prior surgery
surgCors <- cor(bValTopMatrix[,!is.na(as.character(metaData[metaData$prior_surg == 1, "array_sample"]))])
surgCorsTri <- lower.tri(surgCors, diag=TRUE)
surgCors[surgCorsTri] <- NA
surgCors <- as.numeric(surgCors)
surgCors <- surgCors[!is.na(surgCors)]

# put together long table
pairsCors <- data.frame(type = rep("pairs_surgery",3), 
                        cor = c(cor(bValTopMatrix[,"S1_1"],bValTopMatrix[,"S1_2"]),
                                cor(bValTopMatrix[,"S61_1"],bValTopMatrix[,"S61_2"]),
                                cor(bValTopMatrix[,"S3_1"],bValTopMatrix[,"S3_2"])))

pairsCors <- rbind(pairsCors,
                   data.frame(type = rep("pairs_RT",4), 
                              cor = c(cor(bValTopMatrix[,"S6_1"],bValTopMatrix[,"S6_2"]),
                                      cor(bValTopMatrix[,"S9_1"],bValTopMatrix[,"S9_2"]),
                                      cor(bValTopMatrix[,"S60_1"],bValTopMatrix[,"S60_2"]),
                                      cor(bValTopMatrix[,"S61_2"],bValTopMatrix[,"S61_3"]))))

pairsCors <- rbind(pairsCors, data.frame(type = rep("primary",length(priCors)), cor = priCors),
                              data.frame(type = rep("prior_RT",length(rtCors)), cor = rtCors),
                              data.frame(type = rep("prior_surgery",length(surgCors)), cor = surgCors))

pairsCors$type <- factor(pairsCors$type, levels = c("pairs_surgery","pairs_RT","prior_RT","prior_surgery","primary"))

plot <- ggplot(pairsCors) +
  #geom_density(aes(x=cor, fill = type), alpha = 0.5) +
  geom_violin(aes(x = type, y = cor, fill = type), size = 0.25) +
  geom_jitter(aes(x = type, y = cor), size=0.25, alpha = 0.25, shape = 16, position = position_jitter(width = .1)) +
  scale_y_continuous(expand = c(0.01,0.01), lim=c(-1,1.0)) +
  #scale_y_continuous(breaks=seq(0, 0.3, 0.05), expand = c(0.01,0.01), lim=c(0,0.3)) +
  labs(title = "", x="", y="Pearson R") +
  theme_few() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_rect(colour='black'))
ggsave("methylation/pairs/pearson_pairs_violin.pdf", plot, width = 2.75, height = 2, useDingbats = FALSE, limitsize=FALSE)

plot <- ggplot(pairsCors) +
  #geom_density(aes(x=cor, fill = type), alpha = 0.5) +
  geom_boxplot(aes(x = type, y = cor, fill = type), size = 0.25, outlier.size = 0.25) +
  #geom_jitter(aes(x = type, y = cor), size=0.25, alpha = 0.25, shape = 16, position = position_jitter(width = .1)) +
  scale_y_continuous(expand = c(0.01,0.01), lim=c(-1,1.0)) +
  #scale_y_continuous(breaks=seq(0, 0.3, 0.05), expand = c(0.01,0.01), lim=c(0,0.3)) +
  labs(title = "", x="", y="Pearson R") +
  theme_few() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_rect(colour='black'))
ggsave("methylation/pairs/pearson_pairs_boxplot_nojitter.pdf", plot, width = 2.75, height = 2, useDingbats = FALSE, limitsize=FALSE)

######################################################
# Differentially methylated probes between C1 and C2 #
methylclusters <- read.table('/media/liuj/data5/raleigh/vestibular/analysis/methylation/hclust/bval_top2000_v3_cluster_assign_v2.txt',sep='\t', header=FALSE) 
row.names(methylclusters) <- as.character(methylclusters$V1)

# add methylclusters to metadata 
metaData$methylclusters <- methylclusters[row.names(metaData),"V2"]

# get meta for treatment 
metaDmp$methylclusters <- methylclusters[as.character(metaDmp$array_sample),"V2"]

metaDmpTx <- data.frame(metaDmp,stringsAsFactors = FALSE)
metaDmpTx$methylclusters <- as.character(metaDmpTx$methylclusters)
metaDmpTx$methylclusters[metaDmpTx$methylclusters == "c2"] = 1
metaDmpTx$methylclusters[metaDmpTx$methylclusters == "c1"] = 0

dmpTx <- dmpFinder(bVals[,row.names(metaDmpTx)], pheno = metaDmpTx[,"methylclusters"], type="continuous")

# plot fold changes
plot <- ggplot() +
  geom_point(data = dmpTx[dmpTx$qval >= 0.05,], aes(x = beta,y = -log(qval,10)),shape=16,size=0.5,alpha=0.5,color='gray50') + 
  geom_point(data = dmpTx[dmpTx$qval < 0.05,], aes(x = beta,y = -log(qval,10)),shape=16,size=0.5,alpha=0.5,color='green4') + 
  geom_vline(aes(xintercept=0),linetype='dashed',size=0.25) +
  labs(title = paste("C2 vs C1  \n", nrow(dmpSec[dmpTx$qval < 0.05,]),"- probes FDR < 0.05"), x="log odds C2:C1", y="-log10(FDR)") +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_rect(colour='black'))
ggsave("methylation/dmp/plots/C2_vs_C1_05.png", plot, width = 2.2, height = 2.5, units='in', device='png')

# output custom heatmap - probes of interest with gene level info (name and probe coordinate in the rows )
dmpTxSigPos <- dmpTx[dmpTx$beta > 0 & dmpTx$qval < 0.05,]
dmpTxSigNeg <- dmpTx[dmpTx$beta < 0 & dmpTx$qval < 0.05,]
dmpTxSigPos <- dmpTxSigPos[order(-dmpTxSigPos$beta),]
dmpTxSigNeg <- dmpTxSigNeg[order(dmpTxSigNeg$beta),]

probesOI <- c(row.names(dmpTxSigPos[1:1000,]),row.names(dmpTxSigNeg[1:1000,])) # top and bottom 1000 probes by DMP
annEpicSub = annEpic[match(probesOI,annEpic$Name),c(1:4,12:19,24:ncol(annEpic))]
annSimp <- data.frame(annEpicSub[,c(1,2,24)])

# assemble annotated matrix
bValSubMatrix <- bVals[probesOI,]
colnames(bValSubMatrix) <- as.character(targetsSmp[colnames(bValSubMatrix),"sample"])
filler = matrix(data = NA,ncol = 3,nrow = 15)
colnames(filler) <- colnames(annSimp)
row.names(filler) <- c("array_sample","sex","age","size","EOR","preop_growth_rate_percentyr","postop_growth_rate_percentyr","mnp_classifier","mnp_score","methylclusters","prior_rt","prior_surg","clinically_aggressive","LF","LFFP")

bValSubMatrixAnnot <- data.frame(rbind(filler, annSimp[probesOI,]),
                                 rbind(t(metaData[colnames(bValSubMatrix),c("array_sample","sex","age","size","EOR","preop_growth_rate_percentyr","postop_growth_rate_percentyr","mnp_classifier","mnp_score","methylclusters","prior_rt","prior_surg","clinically_aggressive","LF","LFFP")]),
                                       bValSubMatrix))
# output 
write.table(bValSubMatrixAnnot,'methylation/dmp/vs_bValC2C1_dmp_topbot2000_annot.txt', sep = '\t', quote=FALSE, row.names=TRUE, col.names=NA)

# output matrix of ALL significant DMPs with probe annotation 12/28/2018
probesOI <- c(row.names(dmpTxSigPos),row.names(dmpTxSigNeg)) # top and bottom 1000 probes by DMP
annEpicSub = annEpic[match(probesOI,annEpic$Name),c(1:4,12:19,24:ncol(annEpic))]
annSimp <- data.frame(annEpicSub[,c(1,2,24)])

# assemble annotated matrix
bValSubMatrix <- bVals[probesOI,]
colnames(bValSubMatrix) <- as.character(targetsSmp[colnames(bValSubMatrix),"sample"])
filler = matrix(data = NA,ncol = 3,nrow = 15)
colnames(filler) <- colnames(annSimp)
row.names(filler) <- c("array_sample","sex","age","size","EOR","preop_growth_rate_percentyr","postop_growth_rate_percentyr","mnp_classifier","mnp_score","methylclusters","prior_rt","prior_surg","clinically_aggressive","LF","LFFP")

bValSubMatrixAnnot <- data.frame(rbind(filler, annSimp[probesOI,]),
                                 rbind(t(metaData[colnames(bValSubMatrix),c("array_sample","sex","age","size","EOR","preop_growth_rate_percentyr","postop_growth_rate_percentyr","mnp_classifier","mnp_score","methylclusters","prior_rt","prior_surg","clinically_aggressive","LF","LFFP")]),
                                       bValSubMatrix))
# output 
write.table(bValSubMatrixAnnot,'methylation/dmp/vs_bValC2C1_dmp_all_annot.txt', sep = '\t', quote=FALSE, row.names=TRUE, col.names=NA)

# significant DMPs of Raleigh candidate genes 
bValSubMatrixAnnotRaleigh <- bValSubMatrixAnnot[grep(paste(raleighGenes,collapse="|"),bValSubMatrixAnnot$GencodeCompV12_NAME),]
write.table(bValSubMatrixAnnotRaleigh,'methylation/dmp/vs_bValC2C1_dmp_raleigh_annot.txt', sep = '\t', quote=FALSE, row.names=TRUE, col.names=NA)


# correlation of methylation data with RNA-seq
load("/media/liuj/data5/raleigh/vestibular/rnaseq/analysis/vs_rna_kallisto_tpm.rda")

probeKey <- annEpic[,c("chr","pos","GencodeCompV12_NAME")]
probeKey$gene <- sapply(strsplit(probeKey[,3],split=';'),"[",1)

# output probekey in bedformat
probeKeyBed <- data.frame(probeKey)
probeKeyBed <- data.frame(probeKeyBed$chr,probeKeyBed$pos, as.numeric(probeKeyBed$pos + 1), paste(row.names(probeKeyBed),probeKeyBed$gene,sep='-'))
write.table(probeKeyBed,'/home/liuj/genomes/Homo_sapiens/methylation/illuminaEPIC_hg19.bed',row.names=FALSE,col.names=FALSE, sep='\t',quote=FALSE)


for (i in 1:ncol(tpmTblSimpCollapse)){
  smpid = colnames(tpmTblSimpCollapse)[i]
  
  methylvrna <- data.frame(gene = probeKey[row.names(bValsVariance[1:2000,]),"gene"],
                           beta = bVals[row.names(bValsVariance[1:2000,]),as.character(metaData[smpid,"barcode"])],
                           rna = tpmTblSimpCollapse[as.character(probeKey[row.names(bValsVariance[1:2000,]),"gene"]),smpid])
  methylvrna <- methylvrna[!is.na(methylvrna$rna),]
  
  corval = cor(methylvrna$beta,log(methylvrna$rna+1,2))
  
  pdf(paste("methylation/methylvrna/",smpid,"_methylvrna_smooth.pdf",sep=''),width = 2.5,height = 3)
  smoothScatter(methylvrna$beta,log(methylvrna$rna+1,2),xlab = "beta",ylab = "log2(TPM+1)",main = paste(smpid,"R = ",round(corval,3)))
  dev.off()
}

# raleigh genes for methylation
raleighGenes <- as.character(read.table('/media/liuj/data5/raleigh/vestibular/rnaseq/analysis/raleigh_vs_genes.txt',header = FALSE)$V1)
raleighbVals <- data.frame(gene = probeKey[row.names(bVals),"gene"],bVals)
colnames(raleighbVals) <- c("gene",phenodata[colnames(bVals),"sample"])

raleighbVals <- raleighbVals[grep(paste(raleighGenes,collapse="$|^"),raleighbVals$gene),]
row.names(raleighbVals) <- paste(raleighbVals$gene,row.names(raleighbVals),sep='-')
raleighbVals <- raleighbVals[,-1]

raleighbValsAnnot <- rbind(t(metaData[colnames(raleighbVals),c("array_sample","sex","age","size","EOR","preop_growth_rate_percentyr","postop_growth_rate_percentyr","mnp_classifier","mnp_score","prior_rt","prior_surg","clinically_aggressive","LF","LFFP")]),
                         methylClust = t(methylclusters[colnames(raleighbVals),])["V2",],
                         as.matrix(raleighbVals)) # this was log2(bVal+1) previously but turns out it doesn't change much 

write.table(raleighbValsAnnot,'methylation/vs_dna_methylation_raleigh_annot.txt', sep = '\t', quote=FALSE, row.names=TRUE, col.names=NA)

# differentially methylated regions 
phenodmr <- phenodata
phenodmr$methylclusters <- methylclusters[as.character(phenodmr$sample),"V2"]
designMatrix <- model.matrix(~ phenodmr$methylclusters)
dmrs <- bumphunter(MSet.funnorm, design = designMatrix, cutoff = 0.2, B=1000, type="Beta")
dmrs2 <- bumphunter(MSet.funnorm, design = designMatrix, cutoff = 0.4, B=1000, type="Beta") # higher cutoff 

dmrAll <- dmrs$table
dmrSig <- dmrAll[abs(dmrAll$area) > 1,]

# volcano plot
# plot fold changes
plot <- ggplot() +
  geom_point(data = dmrAll[dmrAll$area <= 1,], aes(x = value,y = -log(p.value+(10^-5),10)),shape=16,size=0.5,alpha=0.5,color='gray50') + 
  geom_point(data = dmrAll[dmrAll$area > 1,], aes(x = value,y = -log(p.value+(10^-5),10)),shape=16,size=0.5,alpha=0.5,color='green4') + 
  geom_vline(aes(xintercept=0),linetype='dashed',size=0.25) +
  labs(title = paste("C2 vs C1  \n", nrow(dmrAll[abs(dmrAll$area >= 1),]),"- area > 1"), x="Avg Difference C2:C1", y="-log10(pval)") +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_rect(colour='black'))
ggsave("methylation/dmr/C2_vs_C1_area1.png", plot, width = 2.2, height = 2.5, units='in', device='png')

# output significant DMRs for exploration
write.table(dmrSig,'methylation/dmr/dmrSig_c2c1_area1.txt', sep = '\t', quote=FALSE, row.names=FALSE, col.names=TRUE)

# read back DMR genes comparison to RNAseq 
dmrSigAnnot <- read.table('methylation/dmr/dmrSig_c2c1_area1_closestGenes.bed',sep='\t', header=FALSE) 
dmrSigRNA <- dmrSigAnnot[,c(1:4,8)]

# collapse RNAseq by samples
tpmC1C2 <- data.frame(row.names = row.names(tpmTblSimpCollapse),
                      c2 = apply(tpmTblSimpCollapse[,intersect(colnames(tpmTblSimpCollapse),methylclusters[methylclusters$V2 == "c2","V1"])],1,median),
                      c1 = apply(tpmTblSimpCollapse[,intersect(colnames(tpmTblSimpCollapse),methylclusters[methylclusters$V2 == "c1","V1"])],1,median))
tpmC1C2$c2c1 <- (tpmC1C2$c2+0.1)/(tpmC1C2$c1+0.1)

dmrSigRNA$rnac2c1 <- tpmC1C2[as.character(dmrSigRNA$V8),"c2c1"]
dmrSigRNA <- dmrSigRNA[!is.na(dmrSigRNA$rnac2c1),]

corval = cor(dmrSigRNA$V4, log(dmrSigRNA$rnac2c1,2))
pdf("methylation/dmr/DMRsigvrna_smooth.pdf",width = 2.5,height = 3)
smoothScatter(dmrSigRNA$V4,log(dmrSigRNA$rnac2c1,2),xlab = "Avg C2:C1 DMR difference",ylab = "log2(C2/C1 RNA)",main = paste("R = ",round(corval,3)))
dev.off()
