## Methylation classifier implementation
library(minfi)
library(data.table)
library(xlsx)
library(reshape)
library(caret)
library(dplyr)
library(plyr)
library(kernlab)

setwd('/media/liuj/data5/raleigh/vestibular/analysis')

## Preprocessing of data #################
##########################################
# process data - data set to be analyzed #
baseDir <- "/media/liuj/data5/raleigh/vestibular/DNA_methylation/cytof_runs"

# load array data
targets <- read.metharray.sheet(baseDir)
RGSet <- read.metharray.exp(targets=targets)
RGSet@annotation = c(array= "IlluminaHumanMethylationEPIC", annotation = "ilm10b2.hg19")
phenodata <- pData(RGSet)

# normalize array data
manifest <- getManifest(RGSet)
MSet.funnorm <- preprocessFunnorm(RGSet)

# filter data based on common SNPs, X/Y chromosome, cross reactive probes, detection p-value
# remove SNPs
mSetSqFlt = dropLociWithSnps(MSet.funnorm)

#remove X/Y chromsome
annEpic = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
keep = !(featureNames(mSetSqFlt) %in% annEpic$Name[annEpic$chr %in% c("chrX","chrY")])
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

##########################################
# Apply the classifier 
load("svm_linear_classifier.Rds")

testLines <- t(bVals)
rownames(testLines) <- as.character(phenodata[rownames(testLines),"sample"])
p <- predict(svm_Linear, newdata=testLines)

#######################################
# heatmap comparison to current data ##
datc2c1Clin <- read.table('methylation/dmp/vs_bValC2C1_dmp_topbot2000_annot.txt',sep='\t',header=TRUE)
datc2c1Clin <- datc2c1Clin[16:nrow(datc2c1Clin),c(1,5:ncol(datc2c1Clin))]
row.names(datc2c1Clin) <- as.character(datc2c1Clin$X)
datc2c1Clin <- datc2c1Clin[,-1]

# combine the new methylation data with the older clinical data
dfAdd <- bVals
colnames(dfAdd) <- as.character(phenodata[colnames(dfAdd),"sample"])
datc2c1ClinAdd <- data.frame(datc2c1Clin, dfAdd[row.names(datc2c1Clin),])
datc2c1ClinAdd <- data.frame(row.names = row.names(datc2c1ClinAdd), lapply(datc2c1ClinAdd, as.character), stringsAsFactors=FALSE)

# get metadata back
metaData <- read.xlsx("/home/liuj/Dropbox/Lim_Lab/research_projects/raleigh/vestibular_schwannoma/VS_master_liu.xlsx",1)
metaData <- metaData[!is.na(metaData$array_sample),]
row.names(metaData) <- as.character(metaData$array_sample)
methylclusters <- read.table('/media/liuj/data5/raleigh/vestibular/analysis/methylation/hclust/bval_top2000_v3_cluster_assign_v2.txt',sep='\t', header=FALSE) 
row.names(methylclusters) <- as.character(methylclusters$V1)
metaData$methylclusters <- methylclusters[row.names(metaData),"V2"]
datac2c1ClinMeta <- read.table('methylation/dmp/vs_bValC2C1_dmp_topbot2000_annot.txt',sep='\t',header=TRUE)
datac2c1ClinMeta <- datac2c1ClinMeta[-c(1:15),1:4]
row.names(datac2c1ClinMeta) <- as.character(datac2c1ClinMeta$X)

# add back the metadata that we want
filler = matrix(data = NA,ncol = 4,nrow = 15)
colnames(filler) <- c("X","chr", "pos", "GencodeCompV12_NAME")
row.names(filler) <- c("array_sample","sex","age","size","EOR","preop_growth_rate_percentyr","postop_growth_rate_percentyr","mnp_classifier","mnp_score","methylclusters","prior_rt","prior_surg","clinically_aggressive","LF","LFFP")

# header filler
fillerH <- t(metaData[colnames(datc2c1ClinAdd),c("array_sample","sex","age","size","EOR","preop_growth_rate_percentyr","postop_growth_rate_percentyr","mnp_classifier","mnp_score","methylclusters","prior_rt","prior_surg","clinically_aggressive","LF","LFFP")])
colnames(fillerH) <- c(colnames(datc2c1Clin),colnames(dfAdd))
colnames(fillerH) <- gsub("-",".",colnames(fillerH))

datc2c1ClinAddAnnot <- data.frame(rbind(filler, datac2c1ClinMeta[row.names(datc2c1ClinAdd),]),
                                  rbind(fillerH, datc2c1ClinAdd))
# output 
write.table(datc2c1ClinAddAnnot,'/home/liuj/Dropbox/Lim_Lab/research_projects/raleigh/vestibular_schwannoma/methylation/methyl_cytof/vs_bValC2C1_dmp_topbot2000_addCyTOF.txt', sep = '\t', quote=FALSE, row.names=TRUE, col.names=NA)
