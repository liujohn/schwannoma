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
