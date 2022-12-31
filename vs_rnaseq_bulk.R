# Analysis of RNA-seq bulk data of vestibular schwannoma
library(ggplot2)
library(ggthemes)
library(gplots)
library(GGally)
library(plyr)
library(mixOmics)
library(RColorBrewer)
library(GenomicRanges)
library(data.table)
library(pheatmap)
library(utils)
library(biomaRt)
library(tidyr)
library(xlsx)
library(Seurat)
library(RColorBrewer)
library(DESeq2)
library(WriteXLS)

setwd('/raleighlab/data1/liuj/schwannoma/rnaseq/analysis/')

# prepare data
smpid <- dir(file.path("/media/liuj/data5/raleigh/vestibular/rnaseq/kallisto_out/"))
kaldirs <- file.path("/media/liuj/data5/raleigh/vestibular/rnaseq/kallisto_out", smpid)

# load metadata
metaData <- read.xlsx("/media/liuj/data5/raleigh/vestibular/analysis/VS_master_liu.xlsx", 1)
methylclusters <- read.table('/media/liuj/data5/raleigh/vestibular/analysis/methylation/hclust/bval_top2000_v3_cluster_assign.txt',sep='\t', header=FALSE) 
row.names(methylclusters) <- as.character(methylclusters$V1)
row.names(metaData) <- as.character(metaData$array_sample)

# manual reads of TPM
tpmList <-list()
for (i in 1:length(smpid)){
  df = read.table(paste(kaldirs[[i]],"/abundance.tsv",sep=''),sep='\t',header=TRUE,row.names=1)
  tpmList[[i]] <- df
}
names(tpmList) <- smpid

# TPMs only
tpmTbl <- data.frame(row.names=row.names(tpmList[[1]]), matrix(0, nrow = nrow(tpmList[[1]]), ncol = length(smpid)))
names(tpmTbl) <- smpid

for (i in 1:length(smpid)){
  tpm = tpmList[[i]][row.names(tpmTbl),"tpm"]
  tpmTbl[,i] <- tpm
}

## Make colnames consistent with each other 
colnames(tpmTbl) <- gsub("S33$","S37",colnames(tpmTbl))
colnames(tpmTbl) <- gsub("S1$","S1_1",colnames(tpmTbl))
colnames(tpmTbl) <- gsub("S25$","S36",colnames(tpmTbl))
colnames(tpmTbl) <- gsub("S6$","S6_1",colnames(tpmTbl))
colnames(tpmTbl) <- gsub("S7$","S6_2",colnames(tpmTbl))
colnames(tpmTbl) <- gsub("S9$","S9_2",colnames(tpmTbl))
colnames(tpmTbl) <- gsub("S4$","S3_2",colnames(tpmTbl))
colnames(tpmTbl) <- gsub("S32$","S38",colnames(tpmTbl))

# save TPM table
write.table(tpmTbl,'/media/liuj/data5/raleigh/vestibular/rnaseq/analysis/vs_rna_bulk_TPM.txt', sep = '\t', quote=FALSE, row.names=TRUE, col.names=TRUE)

# collapse TPM matrix by maximum of isoform expression of median across all samples
geneTbl <- read.table(text = row.names(tpmTbl), sep = "|", colClasses = "character") # lengthy run
tpmTblSimp <- data.frame(genes = geneTbl$V6, tpmTbl, median = apply(tpmTbl,1,median))
tpmTblSimpCollapse <- as.data.table(tpmTblSimp)[, .SD[which.max(median)], by=genes]
tpmTblSimpCollapse <- data.frame(row.names=as.character(tpmTblSimpCollapse$genes), tpmTblSimpCollapse)
tpmTblSimpCollapse <- tpmTblSimpCollapse[,2:(ncol(tpmTblSimp)-1)]

# output transformed TPM table
write.table(tpmTblSimpCollapse,'vs_rna_bulk_genes_TPM.txt', sep = '\t', quote=FALSE, row.names=TRUE, col.names=NA)
write.table(round(log(tpmTblSimpCollapse+1,2),4),'vs_rna_bulk_plog2TPM_genes.txt', sep = '\t', quote=FALSE, row.names=TRUE, col.names=NA)

# Get most variable genes using dispersion (variance to mean ratio)
meanDisp <- data.frame(row.names=row.names(tpmTblSimpCollapse), mean = rep(0,nrow(tpmTblSimpCollapse)), variance = rep(0,nrow(tpmTblSimpCollapse)), dispersion = rep(0,nrow(tpmTblSimpCollapse)))
meanDisp$mean <- apply(tpmTblSimpCollapse,1,mean)
meanDisp$variance <- apply(tpmTblSimpCollapse,1,var)
meanDisp$dispersion <- meanDisp$variance/meanDisp$mean
meanDisp <- meanDisp[order(-meanDisp$dispersion),]
plot(log(meanDisp$mean+1),log(meanDisp$dispersion+1),cex=0.5)

# perform PCA
# PCA
pca <- prcomp(t(log(tpmTblSimpCollapse+1,2)),center=TRUE,scale.=FALSE)
prop<-summary(pca)
prop<-prop$importance[2,]

topgenes<-data.frame(pca$rotation[,1:3],maxabs=apply(pca$rotation[,1:3],1,function(x) max(abs(x))))
topgenes<-topgenes[order(-topgenes$maxabs),]
topgenes<-data.frame(topgenes,maxpc=apply(topgenes[,1:3],1,function(x) which.max(abs(x))))
topgenes <- topgenes[topgenes$maxabs > 0,]

# Get PC genes of interest
goi <- row.names((topgenes[order(-topgenes$PC2),]))[1:200]
write.table(goi,'plots/pca/vs_rna_bulk_PC2_top200.txt', sep = '\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

# PCA plotting
smp.pc<-data.frame(pca$x[,1:3], smp=row.names(pca$x), methylclust = as.character(methylclusters[row.names(pca$x),"V2"]))

# colors
scale<-brewer.pal(3,"Dark2")

plot1<-ggplot(smp.pc)+
  geom_point(aes(x=PC1,y=PC2,color=methylclust),size=1)+
  geom_text(aes(x=PC1,y=PC2,label=as.character(smp)),hjust=0.5,vjust=1,size=1)+
  scale_color_manual(name="",values = scale) +
  guides(colour = guide_legend(override.aes = list(size=2)))+
  theme_few() + 
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.ticks=element_line(colour='black',size=0.25),panel.border=element_rect(colour='black'))+ 
  labs(x=paste('PC1 ',prop[1]*100,'%'),y=paste('PC2 ',prop[2]*100,'%'),title='PCA 1 vs 2')
ggsave(paste('plots/pca/pca_1_2.pdf'),plot1,width=3,height=2,useDingbats=FALSE)

plot1<-ggplot(smp.pc)+
  geom_point(aes(x=PC2,y=PC3,color=methylclust),size=1)+
  geom_text(aes(x=PC2,y=PC3,label=as.character(smp)),hjust=0.5,vjust=1,size=1)+
  scale_color_manual(name="",values = scale) +
  guides(colour = guide_legend(override.aes = list(size=2)))+
  theme_few() + 
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.ticks=element_line(colour='black',size=0.25),panel.border=element_rect(colour='black'))+ 
  labs(x=paste('PC2 ',prop[2]*100,'%'),y=paste('PC3 ',prop[3]*100,'%'),title='PCA 2 vs 3')
ggsave(paste('plots/pca/pca_2_3.pdf'),plot1,width=3,height=2,useDingbats=FALSE)


# assembly of heatmap using top variable genes and metadata 
tpmTblSimpCollapseTop <- tpmTblSimpCollapse[row.names(topgenes[1:2000,]),]
tpmTblSimpCollapseTopAnnot <- rbind(t(metaData[colnames(tpmTblSimpCollapseTop),c("array_sample","sex","age","size","EOR","preop_growth_rate_percentyr","postop_growth_rate_percentyr","mnp_classifier","mnp_score","prior_rt","prior_surg","clinically_aggressive","LF","LFFP")]),
                                    methylClust = t(methylclusters[colnames(tpmTblSimpCollapseTop),])["V2",],
                                    as.matrix(log(tpmTblSimpCollapseTop+1,2)))

write.table(tpmTblSimpCollapseTopAnnot,'vs_rna_bulk_top2000PCA_plogTPM_annot.txt', sep = '\t', quote=FALSE, row.names=TRUE, col.names=NA)

tpmTblSimpCollapseTop <- tpmTblSimpCollapse[row.names(meanDisp[1:2000,]),]
tpmTblSimpCollapseTopAnnot <- rbind(t(metaData[colnames(tpmTblSimpCollapseTop),c("array_sample","sex","age","size","EOR","preop_growth_rate_percentyr","postop_growth_rate_percentyr","mnp_classifier","mnp_score","prior_rt","prior_surg","clinically_aggressive","LF","LFFP")]),
                                    methylClust = t(methylclusters[colnames(tpmTblSimpCollapseTop),])["V2",],
                                    as.matrix(log(tpmTblSimpCollapseTop+1,2)))

write.table(tpmTblSimpCollapseTopAnnot,'vs_rna_bulk_top2000dispersion_plogTPM_annot.txt', sep = '\t', quote=FALSE, row.names=TRUE, col.names=NA)

#############################################################
# Differential expression using DESeq2 on hisat alignments  #
counts <- read.table('/media/liuj/data5/raleigh/vestibular/rnaseq/featurecounts_out/vs_rnaseq_counts.txt',header=TRUE,sep='\t',row.names=1, skip=1)
counts.attr <- counts[,1:5]
counts <- counts[,-(1:5)]
names(counts) <- gsub(".bam","",names(counts))
## Make colnames consistent with each other 
colnames(counts) <- gsub("S33$","S37",colnames(counts))
colnames(counts) <- gsub("S1$","S1_1",colnames(counts))
colnames(counts) <- gsub("S25$","S36",colnames(counts))
colnames(counts) <- gsub("S6$","S6_1",colnames(counts))
colnames(counts) <- gsub("S7$","S6_2",colnames(counts))
colnames(counts) <- gsub("S9$","S9_2",colnames(counts))
colnames(counts) <- gsub("S4$","S3_2",colnames(counts))
colnames(counts) <- gsub("S32$","S38",colnames(counts))

genes.attr <- read.table('/home/liuj/genomes/hisat_index/grch38_snp_tran/annotations/genes.attr_table',header=TRUE,row.names=1,sep='\t')

# formulate comparisons
counts.list <- list()
counts.list[["RT_noRT"]]=counts[,c(intersect(row.names(metaData[metaData$prior_rt == 1,]),colnames(counts)),
                                   intersect(row.names(metaData[is.na(metaData$prior_rt),]),colnames(counts)))]
colnames(counts.list[["RT_noRT"]]) <- c(paste("tx", rep(1:length(intersect(row.names(metaData[metaData$prior_rt == 1,]),colnames(counts))),length(intersect(row.names(metaData[metaData$prior_rt == 1,]),colnames(counts))),times=1)),
                                        paste("ctrl", rep(1:length(intersect(row.names(metaData[is.na(metaData$prior_rt),]),colnames(counts))),length(intersect(row.names(metaData[is.na(metaData$prior_rt),]),colnames(counts))),times=1)))

counts.list[["priRec_priNoRec"]]=counts[,c(intersect(row.names(metaData[is.na(metaData$prior_rt) & is.na(metaData$prior_surg) & metaData$LF == 1,]),colnames(counts)),
                                           intersect(row.names(metaData[is.na(metaData$prior_rt) & is.na(metaData$prior_surg) & metaData$LF == 0,]),colnames(counts)))]
colnames(counts.list[["priRec_priNoRec"]]) <- c(paste("tx", rep(1:length(intersect(row.names(metaData[is.na(metaData$prior_rt) & is.na(metaData$prior_surg) & metaData$LF == 1,]),colnames(counts))),length(intersect(row.names(metaData[is.na(metaData$prior_rt) & is.na(metaData$prior_surg) & metaData$LF == 1,]),colnames(counts))),times=1)),
                                                paste("ctrl", rep(1:length(intersect(row.names(metaData[is.na(metaData$prior_rt) & is.na(metaData$prior_surg) & metaData$LF == 0,]),colnames(counts))),length(intersect(row.names(metaData[is.na(metaData$prior_rt) & is.na(metaData$prior_surg) & metaData$LF == 0,]),colnames(counts))),times=1)))

counts.list[["c3_c12"]]=counts[,c(intersect(row.names(methylclusters[methylclusters$V2 == "c3",]),colnames(counts)),
                                  intersect(row.names(methylclusters[methylclusters$V2 != "c3",]),colnames(counts)))]
colnames(counts.list[["c3_c12"]]) <- c(paste("tx", rep(1:length(intersect(row.names(methylclusters[methylclusters$V2 == "c3",]),colnames(counts))),length(intersect(row.names(methylclusters[methylclusters$V2 == "c3",]),colnames(counts))),times=1)),
                                       paste("ctrl", rep(1:length(intersect(row.names(methylclusters[methylclusters$V2 != "c3",]),colnames(counts))),length(intersect(row.names(methylclusters[methylclusters$V2 != "c3",]),colnames(counts))),times=1)))

counts.list[["c3RT_c3noRT"]]=counts[,c(intersect(intersect(row.names(methylclusters[methylclusters$V2 == "c3",]),colnames(counts)),row.names(metaData[metaData$prior_rt == 1,])),
                                       intersect(intersect(row.names(methylclusters[methylclusters$V2 == "c3",]),colnames(counts)),row.names(metaData[is.na(metaData$prior_rt),])))]
colnames(counts.list[["c3RT_c3noRT"]]) <- c(paste("tx", rep(1:length(intersect(intersect(row.names(methylclusters[methylclusters$V2 == "c3",]),colnames(counts)),row.names(metaData[metaData$prior_rt == 1,]))),length(intersect(intersect(row.names(methylclusters[methylclusters$V2 == "c3",]),colnames(counts)),row.names(metaData[metaData$prior_rt == 1,]))),times=1)),
                                            paste("ctrl", rep(1:length(intersect(intersect(row.names(methylclusters[methylclusters$V2 == "c3",]),colnames(counts)),row.names(metaData[is.na(metaData$prior_rt),]))),length(intersect(intersect(row.names(methylclusters[methylclusters$V2 == "c3",]),colnames(counts)),row.names(metaData[is.na(metaData$prior_rt),]))),times=1)))

# DESeq2 results 
res.list <- list()
res.sig.list <- list()

for (i in 1:length(counts.list)){
  sub <- counts.list[[i]]
  cond <- c(rep("tx",length(grep("tx",colnames(sub)))),
            rep("ctrl",length(grep("ctrl",colnames(sub)))))
  dds <- DESeqDataSetFromMatrix(sub,colData=data.frame(sample=colnames(sub),condition=cond),~ condition)
  dds <- DESeq(dds)
  res <- results(dds)
  res <- res[complete.cases(res),]
  res.sig <- res[res$padj < 0.05,] # significance threshold
  res.sig <- res.sig[order(-res.sig$log2FoldChange),]
  res.sig <- data.frame(gname = as.character(genes.attr[row.names(res.sig),"gene_short_name"]), res.sig)
  res <- data.frame(gname = as.character(genes.attr[row.names(res),"gene_short_name"]), res)
  res.list[[names(counts.list)[[i]]]] <- res
  res.sig.list[[names(counts.list)[[i]]]] <- res.sig
}

#volcano plots
for (i in 1:length(res.list)){
  volcano <- data.frame(res.list[[i]])
  volcano <- volcano[order(volcano$log2FoldChange),]
  volcano <- data.frame(gname = genes.attr[row.names(volcano),"gene_short_name"], volcano)
  
  plot <- ggplot() +
    geom_point(data = volcano[volcano$padj > 0.05,], aes(x = log2FoldChange,y = -log(padj,10)),shape=16,size=0.5,alpha=0.5,color='gray50') + 
    geom_point(data = volcano[volcano$padj < 0.05,], aes(x = log2FoldChange,y = -log(padj,10)),shape=16,size=0.5,alpha=0.5,color='green4') + 
    geom_vline(aes(xintercept=0),linetype='dashed',size=0.25) +
    geom_text(data = head(volcano[volcano$padj < 0.05,],10),aes(x = log2FoldChange,y = -log(padj,10), label = gname),vjust=1.5,size=1.5,color='black')+
    geom_text(data = tail(volcano[volcano$padj < 0.05,],10),aes(x = log2FoldChange,y = -log(padj,10), label = gname),vjust=1.5,size=1.5,color='black')+
    labs(title = paste(names(res.list)[[i]],"\n", nrow(volcano[volcano$padj < 0.05,]),"- genes adj p < 0.05"), x="log2 Fold Change", y="-log10(adj p value)") +
    theme_base() +
    theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
          axis.ticks=element_line(colour='black',size=0.25),panel.border=element_rect(colour='black'))
  ggsave(paste("deseq/",names(res.list)[[i]],"_volcano.pdf",sep=''), plot, width = 2.5, height = 2.5, useDingbats = FALSE, limitsize=FALSE)
  
}

WriteXLS(lapply(res.sig.list,data.frame),"deseq/vs_bulk_deseq_sig.xls", SheetNames = names(res.sig.list))
WriteXLS(lapply(res.list,data.frame),"deseq/vs_bulk_deseq_all.xls", SheetNames = names(res.list))

# TPM tables from counts data 
fpkm<-(t(t(counts)/apply(counts,2,sum))*10^6/(counts.attr$Length/1000))
tpm<-t(t(fpkm)/apply(fpkm,2,sum))*10^6
tpm<-tpm[apply(tpm,1,function(x) any(x>0)),]

# generation of heatmap using DE genes subsetting 
tpmRT_noRT <- tpm[row.names(res.sig.list[["RT_noRT"]]),]
row.names(tpmRT_noRT) <- genes.attr[row.names(res.sig.list[["RT_noRT"]]),"gene_short_name"]
tpmRT_noRTAnnot <- rbind(t(metaData[colnames(tpmRT_noRT),c("array_sample","sex","age","size","EOR","preop_growth_rate_percentyr","postop_growth_rate_percentyr","mnp_classifier","mnp_score","prior_rt","prior_surg","clinically_aggressive","LF","LFFP")]),
                    methylClust = t(methylclusters[colnames(tpmRT_noRT),])["V2",],
                    as.matrix(log(tpmRT_noRT+1,2)))

tpmpriRec_priNoRec <- tpm[row.names(res.sig.list[["priRec_priNoRec"]]),]
row.names(tpmpriRec_priNoRec) <- genes.attr[row.names(res.sig.list[["priRec_priNoRec"]]),"gene_short_name"]
tpmpriRec_priNoRecAnnot <- rbind(t(metaData[colnames(tpmpriRec_priNoRec),c("array_sample","sex","age","size","EOR","preop_growth_rate_percentyr","postop_growth_rate_percentyr","mnp_classifier","mnp_score","prior_rt","prior_surg","clinically_aggressive","LF","LFFP")]),
                         methylClust = t(methylclusters[colnames(tpmpriRec_priNoRec),])["V2",],
                         as.matrix(log(tpmpriRec_priNoRec+1,2)))

tpmc3RT_c3noRTc <- tpm[row.names(res.sig.list[["c3RT_c3noRT"]]),]
row.names(tpmc3RT_c3noRTc) <- genes.attr[row.names(res.sig.list[["c3RT_c3noRT"]]),"gene_short_name"]
tpmc3RT_c3noRTAnnot <- rbind(t(metaData[colnames(tpmc3RT_c3noRTc),c("array_sample","sex","age","size","EOR","preop_growth_rate_percentyr","postop_growth_rate_percentyr","mnp_classifier","mnp_score","prior_rt","prior_surg","clinically_aggressive","LF","LFFP")]),
                                 methylClust = t(methylclusters[colnames(tpmc3RT_c3noRTc),])["V2",],
                                 as.matrix(log(tpmc3RT_c3noRTc+1,2)))

tpmc3_c12 <- tpm[row.names(res.sig.list[["c3_c12"]]),]
row.names(tpmc3_c12) <- genes.attr[row.names(res.sig.list[["c3_c12"]]),"gene_short_name"]
tpmc3_c12Annot <- rbind(t(metaData[colnames(tpmc3_c12),c("array_sample","sex","age","size","EOR","preop_growth_rate_percentyr","postop_growth_rate_percentyr","mnp_classifier","mnp_score","prior_rt","prior_surg","clinically_aggressive","LF","LFFP")]),
                             methylClust = t(methylclusters[colnames(tpmc3_c12),])["V2",],
                             as.matrix(log(tpmc3_c12+1,2)))
# output tables
write.table(tpmRT_noRTAnnot,'deseq/vs_rna_bulk_RT_noRT_plogTPM_annot.txt', sep = '\t', quote=FALSE, row.names=TRUE, col.names=NA)
write.table(tpmpriRec_priNoRecAnnot,'deseq/vs_rna_bulk_priRec_priNoRec_plogTPM_annot.txt', sep = '\t', quote=FALSE, row.names=TRUE, col.names=NA)
write.table(tpmc3RT_c3noRTAnnot,'deseq/vs_rna_bulk_c3RT_c3noRT_plogTPM_annot.txt', sep = '\t', quote=FALSE, row.names=TRUE, col.names=NA)
write.table(tpmc3_c12Annot,'deseq/vs_rna_bulk_c3_c12_plogTPM_annot.txt', sep = '\t', quote=FALSE, row.names=TRUE, col.names=NA)

# David's gene table 
raleighGenes <- as.character(read.table('raleigh_vs_genes.txt',header = FALSE)$V1)
raleighTpm <- tpm
row.names(raleighTpm) <- make.unique(as.character(genes.attr[row.names(raleighTpm),"gene_short_name"]))
raleighTpm <- raleighTpm[intersect(raleighGenes,row.names(raleighTpm)),]

raleighTpmAnnot <- rbind(t(metaData[colnames(raleighTpm),c("array_sample","sex","age","size","EOR","preop_growth_rate_percentyr","postop_growth_rate_percentyr","mnp_classifier","mnp_score","prior_rt","prior_surg","clinically_aggressive","LF","LFFP")]),
                             methylClust = t(methylclusters[colnames(raleighTpm),])["V2",],
                             as.matrix(log(raleighTpm+1,2)))

write.table(raleighTpmAnnot,'deseq/vs_rna_bulk_raleigh_plogTPM_annot.txt', sep = '\t', quote=FALSE, row.names=TRUE, col.names=NA)

save(tpmTblSimpCollapse,file='vs_rna_kallisto_tpm.rda')

#########################################################
# differential expression with new methylation clusters #
methylclusters <- read.table('/media/liuj/data5/raleigh/vestibular/analysis/methylation/hclust/bval_top2000_v3_cluster_assign_v2.txt',sep='\t', header=FALSE) 
row.names(methylclusters) <- as.character(methylclusters$V1)

counts.list[["c2RT_c2noRT"]]=counts[,c(intersect(intersect(row.names(methylclusters[methylclusters$V2 == "c2",]),colnames(counts)),row.names(metaData[metaData$prior_rt == 1,])),
                                       intersect(intersect(row.names(methylclusters[methylclusters$V2 == "c2",]),colnames(counts)),row.names(metaData[is.na(metaData$prior_rt),])))]
colnames(counts.list[["c2RT_c2noRT"]]) <- c(paste("tx", rep(1:length(intersect(intersect(row.names(methylclusters[methylclusters$V2 == "c2",]),colnames(counts)),row.names(metaData[metaData$prior_rt == 1,]))),length(intersect(intersect(row.names(methylclusters[methylclusters$V2 == "c2",]),colnames(counts)),row.names(metaData[metaData$prior_rt == 1,]))),times=1)),
                                            paste("ctrl", rep(1:length(intersect(intersect(row.names(methylclusters[methylclusters$V2 == "c2",]),colnames(counts)),row.names(metaData[is.na(metaData$prior_rt),]))),length(intersect(intersect(row.names(methylclusters[methylclusters$V2 == "c2",]),colnames(counts)),row.names(metaData[is.na(metaData$prior_rt),]))),times=1)))

# perform DE
sub <- counts.list[["c2RT_c2noRT"]]
cond <- c(rep("tx",length(grep("tx",colnames(sub)))),
          rep("ctrl",length(grep("ctrl",colnames(sub)))))
dds <- DESeqDataSetFromMatrix(sub,colData=data.frame(sample=colnames(sub),condition=cond),~ condition)
dds <- DESeq(dds)
res <- results(dds)
res <- res[complete.cases(res),]
res.sig <- res[res$padj < 0.05,] # significance threshold
res.sig <- res.sig[order(-res.sig$log2FoldChange),]
res.sig <- data.frame(gname = as.character(genes.attr[row.names(res.sig),"gene_short_name"]), res.sig)
res <- data.frame(gname = as.character(genes.attr[row.names(res),"gene_short_name"]), res)
res.list[["c2RT_c2noRT"]] <- res
res.sig.list[["c2RT_c2noRT"]] <- res.sig

# volcano plot
volcano <- data.frame(res.list[["c2RT_c2noRT"]])
volcano <- volcano[order(volcano$log2FoldChange),]
volcano <- data.frame(gname = genes.attr[row.names(volcano),"gene_short_name"], volcano)

plot <- ggplot() +
  geom_point(data = volcano[volcano$padj > 0.05,], aes(x = log2FoldChange,y = -log(padj,10)),shape=16,size=0.5,alpha=0.5,color='gray50') + 
  geom_point(data = volcano[volcano$padj < 0.05,], aes(x = log2FoldChange,y = -log(padj,10)),shape=16,size=0.5,alpha=0.5,color='green4') + 
  geom_vline(aes(xintercept=0),linetype='dashed',size=0.25) +
  geom_text(data = head(volcano[volcano$padj < 0.05,],10),aes(x = log2FoldChange,y = -log(padj,10), label = gname),vjust=1.5,size=1.5,color='black')+
  geom_text(data = tail(volcano[volcano$padj < 0.05,],10),aes(x = log2FoldChange,y = -log(padj,10), label = gname),vjust=1.5,size=1.5,color='black')+
  labs(title = paste("c2RT_c2noRT\n", nrow(volcano[volcano$padj < 0.05,]),"- genes adj p < 0.05"), x="log2 Fold Change", y="-log10(adj p value)") +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_rect(colour='black'))
ggsave(paste("deseq/c2RT_c2noRT_volcano.pdf",sep=''), plot, width = 2.5, height = 2.5, useDingbats = FALSE, limitsize=FALSE)

tpmc2RT_c2noRT <- tpm[row.names(res.sig.list[["c2RT_c2noRT"]]),]
row.names(tpmc2RT_c2noRT) <- genes.attr[row.names(res.sig.list[["c2RT_c2noRT"]]),"gene_short_name"]
tpmc2RT_c2noRTAnnot <- rbind(t(metaData[colnames(tpmc2RT_c2noRT),c("array_sample","sex","age","size","EOR","preop_growth_rate_percentyr","postop_growth_rate_percentyr","mnp_classifier","mnp_score","prior_rt","prior_surg","clinically_aggressive","LF","LFFP")]),
                        methylClust = t(methylclusters[colnames(tpmc2RT_c2noRT),])["V2",],
                        as.matrix(log(tpmc2RT_c2noRT+1,2)))
# output tables
write.table(tpmc2RT_c2noRTAnnot,'deseq/vs_rna_bulk_c2RT_c2noRT_plogTPM_annot.txt', sep = '\t', quote=FALSE, row.names=TRUE, col.names=NA)

##############################################
# differential expression using new C2 vs C1 #
counts.list[["c2_c1"]]=counts[,c(intersect(row.names(methylclusters[methylclusters$V2 == "c2",]),colnames(counts)),
                                  intersect(row.names(methylclusters[methylclusters$V2 != "c2",]),colnames(counts)))]
colnames(counts.list[["c2_c1"]]) <- c(paste("tx", rep(1:length(intersect(row.names(methylclusters[methylclusters$V2 == "c2",]),colnames(counts))),length(intersect(row.names(methylclusters[methylclusters$V2 == "c2",]),colnames(counts))),times=1)),
                                       paste("ctrl", rep(1:length(intersect(row.names(methylclusters[methylclusters$V2 != "c2",]),colnames(counts))),length(intersect(row.names(methylclusters[methylclusters$V2 != "c2",]),colnames(counts))),times=1)))

# perform DE
sub <- counts.list[["c2_c1"]]
cond <- c(rep("tx",length(grep("tx",colnames(sub)))),
          rep("ctrl",length(grep("ctrl",colnames(sub)))))
dds <- DESeqDataSetFromMatrix(sub,colData=data.frame(sample=colnames(sub),condition=cond),~ condition)
dds <- DESeq(dds)
res <- results(dds)
res <- res[complete.cases(res),]
res.sig <- res[res$padj < 0.05,] # significance threshold
res.sig <- res.sig[order(-res.sig$log2FoldChange),]
res.sig <- data.frame(gname = as.character(genes.attr[row.names(res.sig),"gene_short_name"]), res.sig)
res <- data.frame(gname = as.character(genes.attr[row.names(res),"gene_short_name"]), res)
res.list[["c2_c1"]] <- res
res.sig.list[["c2_c1"]] <- res.sig

write.table(res.sig.list[["c2_c1"]],'deseq/c2_c1_deseq_res_sig.txt', sep = '\t', quote=FALSE, row.names=TRUE, col.names=NA)


# volcano plot
volcano <- data.frame(res.list[["c2_c1"]])
volcano <- volcano[order(volcano$log2FoldChange),]
volcano <- data.frame(gname = genes.attr[row.names(volcano),"gene_short_name"], volcano)

plot <- ggplot() +
  geom_point(data = volcano[volcano$padj > 0.05,], aes(x = log2FoldChange,y = -log(padj,10)),shape=16,size=0.5,alpha=0.5,color='gray50') + 
  geom_point(data = volcano[volcano$padj < 0.05,], aes(x = log2FoldChange,y = -log(padj,10)),shape=16,size=0.5,alpha=0.5,color='green4') + 
  geom_vline(aes(xintercept=0),linetype='dashed',size=0.25) +
  geom_text(data = head(volcano[volcano$padj < 0.05,],10),aes(x = log2FoldChange,y = -log(padj,10), label = gname),vjust=1.5,size=1.5,color='black')+
  geom_text(data = tail(volcano[volcano$padj < 0.05,],10),aes(x = log2FoldChange,y = -log(padj,10), label = gname),vjust=1.5,size=1.5,color='black')+
  labs(title = paste("c2_c1\n", nrow(volcano[volcano$padj < 0.05,]),"- genes adj p < 0.05"), x="log2 Fold Change", y="-log10(adj p value)") +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_rect(colour='black'))
ggsave(paste("deseq/c2_c1_volcano.pdf",sep=''), plot, width = 2.5, height = 2.5, useDingbats = FALSE, limitsize=FALSE)

tpmc2_c1 <- tpm[row.names(res.sig.list[["c2_c1"]]),]
row.names(tpmc2_c1) <- genes.attr[row.names(res.sig.list[["c2_c1"]]),"gene_short_name"]
tpmc2_c1Annot <- rbind(t(metaData[colnames(tpmc2_c1),c("array_sample","sex","age","size","EOR","preop_growth_rate_percentyr","postop_growth_rate_percentyr","mnp_classifier","mnp_score","prior_rt","prior_surg","clinically_aggressive","LF","LFFP")]),
                             methylClust = t(methylclusters[colnames(tpmc2_c1),])["V2",],
                             as.matrix(log(tpmc2_c1+1,2)))
# output tables
write.table(tpmc2_c1Annot,'deseq/vs_rna_bulk_c2_c1_plogTPM_annot.txt', sep = '\t', quote=FALSE, row.names=TRUE, col.names=NA)

#########################################################
# differential expression C2 no RT vs C1 #
counts.list[["c2noRT_c1"]]=counts[,c(intersect(intersect(row.names(methylclusters[methylclusters$V2 == "c2",]),colnames(counts)),row.names(metaData[is.na(metaData$prior_rt),])),
                                     intersect(row.names(methylclusters[methylclusters$V2 != "c2",]),colnames(counts)))]
colnames(counts.list[["c2noRT_c1"]]) <- c(paste("tx", rep(1:length(intersect(intersect(row.names(methylclusters[methylclusters$V2 == "c2",]),colnames(counts)),row.names(metaData[is.na(metaData$prior_rt),]))),length(intersect(intersect(row.names(methylclusters[methylclusters$V2 == "c2",]),colnames(counts)),row.names(metaData[is.na(metaData$prior_rt),]))),times=1)),
                                            paste("ctrl", rep(1:length(intersect(row.names(methylclusters[methylclusters$V2 != "c2",]),colnames(counts)))),length(intersect(row.names(methylclusters[methylclusters$V2 != "c2",]),colnames(counts))),times=1))

# perform DE
sub <- counts.list[["c2noRT_c1"]]
cond <- c(rep("tx",length(grep("tx",colnames(sub)))),
          rep("ctrl",length(grep("ctrl",colnames(sub)))))
dds <- DESeqDataSetFromMatrix(sub,colData=data.frame(sample=colnames(sub),condition=cond),~ condition)
dds <- DESeq(dds)
res <- results(dds)
res <- res[complete.cases(res),]
res.sig <- res[res$padj < 0.05,] # significance threshold
res.sig <- res.sig[order(-res.sig$log2FoldChange),]
res.sig <- data.frame(gname = as.character(genes.attr[row.names(res.sig),"gene_short_name"]), res.sig)
res <- data.frame(gname = as.character(genes.attr[row.names(res),"gene_short_name"]), res)
res.list[["c2noRT_c1"]] <- res
res.sig.list[["c2noRT_c1"]] <- res.sig

# volcano plot
volcano <- data.frame(res.list[["c2noRT_c1"]])
volcano <- volcano[order(volcano$log2FoldChange),]
volcano <- data.frame(gname = genes.attr[row.names(volcano),"gene_short_name"], volcano)

plot <- ggplot() +
  geom_point(data = volcano[volcano$padj > 0.05,], aes(x = log2FoldChange,y = -log(padj,10)),shape=16,size=0.5,alpha=0.5,color='gray50') + 
  geom_point(data = volcano[volcano$padj < 0.05,], aes(x = log2FoldChange,y = -log(padj,10)),shape=16,size=0.5,alpha=0.5,color='green4') + 
  geom_vline(aes(xintercept=0),linetype='dashed',size=0.25) +
  geom_text(data = head(volcano[volcano$padj < 0.05,],10),aes(x = log2FoldChange,y = -log(padj,10), label = gname),vjust=1.5,size=1.5,color='black')+
  geom_text(data = tail(volcano[volcano$padj < 0.05,],10),aes(x = log2FoldChange,y = -log(padj,10), label = gname),vjust=1.5,size=1.5,color='black')+
  labs(title = paste("c2noRT_c1\n", nrow(volcano[volcano$padj < 0.05,]),"- genes adj p < 0.05"), x="log2 Fold Change", y="-log10(adj p value)") +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_rect(colour='black'))
ggsave(paste("deseq/c2noRT_c1_volcano.pdf",sep=''), plot, width = 2.5, height = 2.5, useDingbats = FALSE, limitsize=FALSE)

tpmc2noRT_c1 <- tpm[row.names(res.sig.list[["c2noRT_c1"]]),]
row.names(tpmc2noRT_c1) <- genes.attr[row.names(res.sig.list[["c2noRT_c1"]]),"gene_short_name"]
tpmc2noRT_c1Annot <- rbind(t(metaData[colnames(tpmc2noRT_c1),c("array_sample","sex","age","size","EOR","preop_growth_rate_percentyr","postop_growth_rate_percentyr","mnp_classifier","mnp_score","prior_rt","prior_surg","clinically_aggressive","LF","LFFP")]),
                             methylClust = t(methylclusters[colnames(tpmc2noRT_c1),])["V2",],
                             as.matrix(log(tpmc2noRT_c1+1,2)))
# output tables
write.table(tpmc2noRT_c1Annot,'deseq/vs_rna_bulk_c2noRT_c1_plogTPM_annot.txt', sep = '\t', quote=FALSE, row.names=TRUE, col.names=NA)

