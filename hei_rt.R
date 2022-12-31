# Analysis of RNA-seq bulk data of HEI193 post RT 
library(ggplot2)
library(ggthemes)
library(GGally)
library(reshape)
library(data.table)
library(DESeq2)
library(WriteXLS)
library(pheatmap)
library(RColorBrewer)

setwd('/raleighlab/data1/liuj/schwannoma/rnaseq_cells/analysis')

## read in data 
counts <- read.table('/media/liuj/data5/raleigh/vestibular/rnaseq_cells/featurecounts_out/HEI_frnaseq_counts.txt',header=TRUE,sep='\t',row.names=1, skip=1)
counts.attr <- counts[,1:5]
counts <- counts[,-(1:5)]
names(counts) <- gsub(".bam.bam","",names(counts))
names(counts) <- gsub("HEI195_1","HEI_0Gy_1",names(counts))
names(counts) <- gsub("HEI195_2","HEI_0Gy_2",names(counts))
names(counts) <- gsub("HEI195_3","HEI_0Gy_3",names(counts))
names(counts) <- gsub("HEI195_4","HEI_12.5Gy_5d_1",names(counts))
names(counts) <- gsub("HEI195_5","HEI_12.5Gy_5d_2",names(counts))
names(counts) <- gsub("HEI195_6","HEI_12.5Gy_5d_3",names(counts))
names(counts) <- gsub("HEI195_7","HEI_12.5Gy_19d_1",names(counts))
names(counts) <- gsub("HEI195_8","HEI_12.5Gy_19d_2",names(counts))
names(counts) <- gsub("HEI195_9","HEI_12.5Gy_19d_3",names(counts))

genes.attr <- read.table('/home/liuj/genomes/hisat_index/grch38_snp_tran/annotations/genes.attr_table',header=TRUE,row.names=1,sep='\t')

# TPM table
fpkm<-(t(t(counts)/apply(counts,2,sum))*10^6/(counts.attr$Length/1000))
tpm<-t(t(fpkm)/apply(fpkm,2,sum))*10^6
tpm<-tpm[apply(tpm,1,function(x) any(x>0)),]

tpmAnnot <- tpm[order(-apply(tpm,1,mean)),]
row.names(tpmAnnot) <- make.unique(as.character(genes.attr[row.names(tpmAnnot),"gene_short_name"]))

save(tpmAnnot,file = "HEI_rnaseq_TPM.rds")

# DESeq normalize
ncounts <- t(t(counts)/estimateSizeFactorsForMatrix(counts))

# some QC ----
counts0.1 <- counts[apply(counts,1,function(x) any(x>0.1)),]
ggpairs(log(counts0.1+1,2))

pdf('plots/hei_cor_heatmap_allgenes.pdf', width=4, height=2)
breaklist=seq(0.9,1,0.001)
pheatmap(cor(log(counts0.1+1,2)),cluster_rows = TRUE, cluster_cols = TRUE, color = colorRampPalette((brewer.pal(n = 7, name = "OrRd")))(length(breaklist)),
         breaks = breaklist,show_colnames = FALSE, border_color=NA,cellwidth=10,cellheight=10,treeheight_row=10,treeheight_col=10,fontsize=8,
         main="HEI193 pearson R - all genes")
dev.off()

# formulate comparisons
counts.list <- list()
counts.list[["HEI_RT_5d"]]=counts[,c(c("HEI_12.5Gy_5d_1", "HEI_12.5Gy_5d_2","HEI_12.5Gy_5d_3"),c("HEI_0Gy_1", "HEI_0Gy_2", "HEI_0Gy_3"))]
counts.list[["HEI_RT_19d"]]=counts[,c(c("HEI_12.5Gy_19d_1", "HEI_12.5Gy_19d_2","HEI_12.5Gy_19d_3"),c("HEI_0Gy_1", "HEI_0Gy_2", "HEI_0Gy_3"))]
counts.list[["HEI_RT_19d_v_5d"]]=counts[,c(c("HEI_12.5Gy_19d_1", "HEI_12.5Gy_19d_2","HEI_12.5Gy_19d_3"),c("HEI_12.5Gy_5d_1", "HEI_12.5Gy_5d_2","HEI_12.5Gy_5d_3"))]

res.list <- list()
res.sig.list <- list()

for (i in 1:length(counts.list)){
  sub <- counts.list[[i]]
  cond <- c(rep("kd",3),
            rep("ctrl",3))
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
    xlim(-10, 10) +
    #geom_text(data = head(volcano[volcano$padj < 0.05,],50),aes(x = log2FoldChange,y = -log(padj,10), label = gname),vjust=1.5,size=1.5,color='black')+
    #geom_text(data = tail(volcano[volcano$padj < 0.05,],50),aes(x = log2FoldChange,y = -log(padj,10), label = gname),vjust=1.5,size=1.5,color='black')+
    labs(title = paste(names(res.list)[[i]],"\n", nrow(volcano[volcano$padj < 0.05,]),"- genes adj p < 0.05"), x="log2 Fold Change", y="-log10(adj p value)") +
    theme_base() +
    theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
          axis.ticks=element_line(colour='black',size=0.25),panel.border=element_rect(colour='black'))
  ggsave(paste("plots/", names(res.list)[[i]],"_volcanonolabel.pdf",sep=''), plot, width = 2.5, height = 2.5, useDingbats = FALSE, limitsize=FALSE)
  
}

WriteXLS(lapply(res.sig.list,data.frame),"hei_rt_deseq_sig.xls", SheetNames = names(res.sig.list))
WriteXLS(lapply(res.list,data.frame),"hei_rt_deseq_all.xls", SheetNames = names(res.list))

# Heatmap comparison to qPCR
genes <- c('EGR2','APOE','APOD','JUN','MPZ','SOX6','SOX10','GLI1','APOC1','EGR1','PTCH1')
goiD5 <- res.list[['HEI_RT_5d']][grep(paste(genes,collapse='$|^'),res.list[['HEI_RT_5d']][,"gname"]),]
goiD19 <- res.list[['HEI_RT_19d']][grep(paste(genes,collapse='$|^'),res.list[['HEI_RT_19d']][,"gname"]),]
row.names(goiD5) <- as.character(goiD5$gname)
row.names(goiD19) <- as.character(goiD19$gname)

fcPanel <- data.frame(gene = genes, d5 = goiD5[genes,"log2FoldChange"], d19 = goiD19[genes,"log2FoldChange"])
row.names(fcPanel) <- as.character(fcPanel$gene)
fcPanel <- fcPanel[,-1]

pdf('plots/hei_genepanel_heatmap.pdf', width=4, height=3)
breaklist=seq(-6,6,0.01)
pheatmap(fcPanel,cluster_rows = FALSE, cluster_cols = FALSE, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breaklist)),
         breaks = breaklist,show_colnames = TRUE, border_color=NA,cellwidth=10,cellheight=10,treeheight_row=10,treeheight_col=10,fontsize=8,
         main="HEI193 Gene Panel")
dev.off()

# Correlations using DE genes
deGenes <- unique(row.names(res.sig[1]),row.names(res.sig[2]))
counts0.1DE <- counts0.1[deGenes,]

pdf('plots/hei_cor_heatmap_DEgenes.pdf', width=4, height=2)
breaklist=seq(0.8,1,0.001)
pheatmap(cor(log(counts0.1DE+1,2)),cluster_rows = TRUE, cluster_cols = TRUE, color = colorRampPalette((brewer.pal(n = 7, name = "OrRd")))(length(breaklist)),
         breaks = breaklist,show_colnames = FALSE, border_color=NA,cellwidth=10,cellheight=10,treeheight_row=10,treeheight_col=10,fontsize=8,
         main="HEI193 pearson R - DE genes")
dev.off()

# gene panel #2
genes <- c('APOE','APOD','APOC1','IL1A','CCL3','SOX6','SOX8','SOX11','MPZL1','GLI2','PTPRG','TNC','FOXM1')
# heatmap of metabolism genes 11/21/22
genes <- c("FLAD1", "IDH1", "IDH2", "IDH3", "MDH1", "MDH2", "GULOP", "GULO", "SDH", "CYP24A1", "CYP27B1")

goiD5 <- res.list[['HEI_RT_5d']][grep(paste(genes,collapse='$|^'),res.list[['HEI_RT_5d']][,"gname"]),]
goiD19 <- res.list[['HEI_RT_19d']][grep(paste(genes,collapse='$|^'),res.list[['HEI_RT_19d']][,"gname"]),]
row.names(goiD5) <- as.character(goiD5$gname)
row.names(goiD19) <- as.character(goiD19$gname)

fcPanel <- data.frame(gene = genes, d5 = goiD5[genes,"log2FoldChange"], d19 = goiD19[genes,"log2FoldChange"])
row.names(fcPanel) <- as.character(fcPanel$gene)
fcPanel <- fcPanel[,-1]

pdf('plots/hei_genepanel_metab_heatmap.pdf', width=4, height=5)
breaklist=seq(-3,3,0.01)
pheatmap(fcPanel,cluster_rows = FALSE, cluster_cols = FALSE, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breaklist)),
         breaks = breaklist,show_colnames = TRUE, border_color=NA,cellwidth=10,cellheight=10,treeheight_row=10,treeheight_col=10,fontsize=8,
         main="HEI193 Gene Panel Metabolism")
dev.off()
write.table(fcPanel,'plots/hei_genepanel_metab_heatmap.txt',sep='\t',quote=FALSE,row.names=TRUE,col.names=NA)


# heatmap of qPCR data
hei_qPCR <- read.table('/home/liuj/Dropbox/Lim_Lab/research_projects/raleigh/vestibular_schwannoma/qPCR/vs_qpcr_12.5Gy_means.txt', sep='\t', row.names=1, header=TRUE)
hsc_qPCR <- read.table('/home/liuj/Dropbox/Lim_Lab/research_projects/raleigh/vestibular_schwannoma/qPCR/hsc_qpcr_12.5Gy_means.txt', sep='\t', row.names=1, header=TRUE)

pdf('/home/liuj/Dropbox/Lim_Lab/research_projects/raleigh/vestibular_schwannoma/qPCR/hsc_qpcr_12.5Gy_heatmap.pdf', width=4, height=5)
breaklist=seq(-2,2,0.01)
pheatmap(hsc_qPCR,cluster_rows = FALSE, cluster_cols = FALSE, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breaklist)),
         breaks = breaklist,show_colnames = TRUE, border_color=NA,cellwidth=10,cellheight=10,treeheight_row=10,treeheight_col=10,fontsize=8,
         main="HSC")
dev.off()

pdf('/home/liuj/Dropbox/Lim_Lab/research_projects/raleigh/vestibular_schwannoma/qPCR/vs_qpcr_12.5Gy_heatmap.pdf', width=4, height=5)
breaklist=seq(-6,6,0.01)
pheatmap(hei_qPCR,cluster_rows = FALSE, cluster_cols = FALSE, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breaklist)),
         breaks = breaklist,show_colnames = TRUE, border_color=NA,cellwidth=10,cellheight=10,treeheight_row=10,treeheight_col=10,fontsize=8,
         main="HEI193")
dev.off()
