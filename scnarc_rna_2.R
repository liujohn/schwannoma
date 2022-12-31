## Analysis of HEI-193 snARC-seq combined - Part 1
library(Signac)
library(Seurat)
library(plyr)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(ggseqlogo)
library(ggforce)
library(biovizBase)
library(qlcMatrix)
library(ggplot2)
library(ggthemes)
library(plyr)
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(future)
library(pheatmap)
library(rGREAT)
library(DESeq2)
library(MAST)
library(cicero)
library(SeuratWrappers)
library(monocle3)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)
library(genomation)

setwd('/raleighlab/data1/liuj/schwannoma/snarc-seq/analysis/snarc_rna_2')

# load all vs in vivo RNA data
smpid <- dir(file.path("/raleighlab/data1/liuj/schwannoma/snarc-seq/cellranger_arc_out/"))

# get gene annotations for hg38
#annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86,standard.chromosomes = TRUE)
#names(annotation@seqinfo) <- paste("chr",names(annotation@seqinfo),sep='')

# read external annotations from Jason 
annot <- readRDS(file = '../HS_ENsg_annotations')

# build seurat objects
dataGene.list <- list()

for (i in 1:length(smpid)){
  df <- Read10X(data.dir = paste("/raleighlab/data1/liuj/schwannoma/snarc-seq/cellranger_arc_out/",smpid[i],"/outs/filtered_feature_bc_matrix/",sep=''))
  dataGene.list[[i]] <- CreateSeuratObject(counts=df$`Gene Expression`, project = smpid[i], assay = "RNA")
  dataGene.list[[i]][["ATAC"]] <- CreateChromatinAssay(counts = df$Peaks, sep = c(":", "-"), 
                                                       fragments = paste("/raleighlab/data1/liuj/schwannoma/snarc-seq/cellranger_arc_out/",smpid[i],"/outs/atac_fragments.tsv.gz",sep=''), 
                                                       annotation = annot)
  dataGene.list[[i]]$origsamplesent <- smpid[i]
}

names(dataGene.list) <- smpid

dat.all <- merge(x = dataGene.list[[1]], y = dataGene.list[2:length(dataGene.list)])

# plot QC metrics
DefaultAssay(dat.all) <- "ATAC"

dat.all <- NucleosomeSignal(dat.all)
dat.all <- TSSEnrichment(dat.all, fast=FALSE, assay = 'ATAC',verbose=TRUE)

VlnPlot(
  object = dat.all,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)
ggsave(filename = "qc_violin.pdf",width = 16, height=10)

# qc subset 
dat.sub <- subset(
  x = dat.all,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nFeature_RNA > 50 &
    nucleosome_signal < 4 &
    TSS.enrichment > 1
)

dat.sub

peaks <- CallPeaks(dat.sub, macs2.path = "/c4/home/liuj/anaconda3/envs/py27/bin/macs2")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(dat.sub),
  features = peaks,
  cells = colnames(dat.sub)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
dat.sub[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = Fragments(dat.sub),
  annotation = annot
)

DefaultAssay(dat.sub) <- "RNA"
dat.sub <- SCTransform(dat.sub)
dat.sub <- RunPCA(dat.sub)

DefaultAssay(dat.sub) <- "peaks"
dat.sub <- FindTopFeatures(dat.sub, min.cutoff = 5)
dat.sub <- RunTFIDF(dat.sub)
dat.sub <- RunSVD(dat.sub)

dat.sub <- FindMultiModalNeighbors(
  object = dat.sub,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
dat.sub <- RunUMAP(
  object = dat.sub,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

DimPlot(dat.sub, label = TRUE, repel = TRUE, reduction = "umap") + NoLegend()
ggsave(filename = "vs_snarc_dimplot.pdf",width = 8, height=6)

DefaultAssay(dat.sub) <- "peaks"

# first compute the GC content for each peak
dat.sub <- RegionStats(dat.sub, genome = BSgenome.Hsapiens.UCSC.hg38)

# link peaks to genes
dat.sub <- LinkPeaks(
  object = dat.sub,
  peak.assay = "peaks",
  expression.assay = "SCT",
  #genes.use = c("LYZ", "MS4A1")
)

goi <- "FGF18"
CoveragePlot(dat.sub, region = goi,features = goi)
ggsave(filename = paste("vs_snarc_coverage_",goi,".pdf", sep=''),width = 8, height=6)

TilePlot(dat.sub, region = goi,tile.cells = 10)

# analysis of linked peaks 
dat.sub.links <- Links(dat.sub)
dat.sub.links <- dat.sub.links[order(-dat.sub.links$zscore),]
head(dat.sub.links)

# load in CROP-seq annotations
guideCallsTbl <- read.table('../ArchR_guide_calls.csv',header=TRUE,sep=',')
guideCallsTbl$suffix <- guideCallsTbl$samples

df <- data.frame(dat.sub$orig.ident, str_split_fixed(colnames(dat.sub), "_", 2)[,2])
df <- unique(df)
row.names(df) <- as.character(df$dat.sub.orig.ident)

# make replacements 
guideCallsTbl$suffix <- gsub("NoRT_1","HEI193_0_1",guideCallsTbl$suffix)
guideCallsTbl$suffix <- gsub("NoRT_2","HEI193_0_2",guideCallsTbl$suffix)
guideCallsTbl$suffix <- gsub("NoRT_3","HEI193_0_3",guideCallsTbl$suffix)
guideCallsTbl$suffix <- gsub("12_Gy_1","HEI193_12_5_1",guideCallsTbl$suffix)
guideCallsTbl$suffix <- gsub("12_Gy_2","HEI193_12_5_2",guideCallsTbl$suffix)
guideCallsTbl$suffix <- gsub("5x1.8_Gy_1","HEI193_5x1_8_1",guideCallsTbl$suffix)
guideCallsTbl$suffix <- gsub("5x1.8_Gy_2","HEI193_5x1_8_2",guideCallsTbl$suffix)
guideCallsTbl$suffix <- gsub("5x1.8_Gy_3","HEI193_5x1_8_3",guideCallsTbl$suffix)

guideCallsTbl$suffix2 <- paste("-1_",df[as.character(guideCallsTbl$suffix),2],sep='')
row.names(guideCallsTbl) <- paste(guideCallsTbl$barcode,guideCallsTbl$suffix2,sep='')

guideCallsTbl$gene <- str_split_fixed(guideCallsTbl$guide, "_", 3)[,1]
guideCallsTbl$geneCond <- paste(guideCallsTbl$gene,guideCallsTbl$samples,sep='_')
guideCallsTbl$geneCondShort <- paste(str_split_fixed(guideCallsTbl$geneCond, "_", 3)[,1],
                                     str_split_fixed(guideCallsTbl$geneCond, "_", 3)[,2],sep="_")

guideCallsTbl$cond <- str_split_fixed(guideCallsTbl$geneCondShort, "_", 2)[,2]

# get summary data - only those also with ATAC data 
df <- guideCallsTbl[intersect(colnames(dat.sub),colnames(dat.sub)),]
df <- melt(table(df[,c("gene","cond")]))

plot <- ggplot(df, aes(x = gene,cond)) +
  geom_tile(aes(fill = value)) +
  geom_text(aes(label = value),size=1) +
  scale_fill_gradient(low = "white", high = "steelblue",limit=c(0,100),na.value="steelblue") +
  labs(title = "Number of cells with sgRNAs detected", x="", y="") +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_blank(),panel.background=element_blank(),axis.line = element_line(size=0.25))
ggsave("vs_snarc_sgRNAs_heatmap_sub.png", plot, width =5.5, height = 1.5)

# combine sgRNA annotations with main data 
dat.sub$sgRNACond <- guideCallsTbl[colnames(dat.sub),"geneCondShort"]

# differential peaks within conditions WITH downsampling 
Idents(dat.sub) <- dat.sub$sgRNACond
DefaultAssay(dat.sub) <- 'peaks'

resList <- list()
smpid <- unique(dat.sub$sgRNACond)[grep("NoRT",unique(dat.sub$sgRNACond))]
smpid <- smpid[grep("non-targeting",invert=TRUE,smpid)]

for (i in 1:length(smpid)){
    dnsmp = ncol(subset(dat.sub, idents=smpid[i])) # downsample cells
    #res <- FindMarkers(object = dat.sub,ident.1 = "non-targeting_NoRT", ident.2 = smpid[i],min.pct = 0.05,test.use = 'LR',max.cells.per.ident = dnsmp)
    #res <- FindMarkers(object = dat.sub,ident.1 = "non-targeting_NoRT", ident.2 = smpid[i],min.pct = 0.05,test.use = 'wilcox',max.cells.per.ident = dnsmp)
    res <- FindMarkers(object = dat.sub,ident.1 = "non-targeting_NoRT", ident.2 = smpid[i],min.pct = 0.05,test.use = 'MAST',max.cells.per.ident = dnsmp)
    resList[[i]] <- res 
    genes.to.label = c(row.names(head(res,30)),row.names(tail(res,30)))
    
    #plot it 
    p1 <- ggplot(res, aes(avg_log2FC,-log(p_val,10))) + geom_point(color='red2') + ggtitle(smpid[i]) +geom_vline(xintercept = 0) + theme(legend.position = "none")
    p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE,xnudge=0,ynudge=0)
    ggsave(paste("snARC_RT_noRT_wilcox",smpid[i],"_scatter.png",sep=''), p1, width = 6, height = 6)
  }

names(resList) <- smpid[1:55] # last few had errors

# now do for RT 
resListRT <- list()
smpid <- unique(dat.sub$sgRNACond)[grep("5x1.8",unique(dat.sub$sgRNACond))]
smpid <- smpid[grep("non-targeting",invert=TRUE,smpid)]

for (i in 1:length(smpid)){
  dnsmp = ncol(subset(dat.sub, idents=smpid[i])) # downsample cells
  #res <- FindMarkers(object = dat.sub,ident.1 = "non-targeting_5x1.8", ident.2 = smpid[i],min.pct = 0.05,test.use = 'LR',max.cells.per.ident = dnsmp)
  #res <- FindMarkers(object = dat.sub,ident.1 = "non-targeting_5x1.8", ident.2 = smpid[i],min.pct = 0.05,test.use = 'wilcox',max.cells.per.ident = dnsmp)
  res <- FindMarkers(object = dat.sub,ident.1 = "non-targeting_5x1.8", ident.2 = smpid[i],min.pct = 0.05,test.use = 'MAST',max.cells.per.ident = dnsmp)
  resListRT[[i]] <- res 
  genes.to.label = c(row.names(head(res,30)),row.names(tail(res,30)))
  
  #plot it 
  p1 <- ggplot(res, aes(avg_log2FC,-log(p_val,10))) + geom_point(color='red2') + ggtitle(smpid[i]) +geom_vline(xintercept = 0) + theme(legend.position = "none")
  p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE,xnudge=0,ynudge=0)
  ggsave(paste("snARC_RT_noRT__wilcox",smpid[i],"_scatter.png",sep=''), p1, width = 6, height = 6)
}

names(resListRT) <- smpid

resListAll <- c(resList,resListRT)

resListAll.MAST <- resListAll
#resListAll.wilcox <- resListAll
#resListAll.LR <- resListAll # backup if i run DARs with wilcox

# DARS between noRT and RT NTC
res <- FindMarkers(object = dat.sub,ident.1 = "non-targeting_NoRT", ident.2 = "non-targeting_5x1.8",min.pct = 0.05,test.use = 'LR')
res$gene <- annot[nearest(StringToGRanges(row.names(res)), annot, ignore.strand=TRUE)]$gene_name
res <- res[order(res$avg_log2FC),]
row.names(res) <- make.unique(as.character(res$gene))

#plot it 
genes.to.label = c(row.names(head(res,30)),row.names(tail(res,30)))
p1 <- ggplot(res, aes(avg_log2FC,-log(p_val,10))) + geom_point(color='red2') + ggtitle("RT vs no RT NTC") +geom_vline(xintercept = 0) + theme(legend.position = "none")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE,xnudge=0,ynudge=0,max.overlaps = 20)
ggsave("snARC_RT_vs_noRT_NTC_volcano.png", p1, width = 6, height = 6)

# remove 12.5 Gy samples
ids <- as.character(unique((Idents(dat.sub))))[grep("12$",as.character(unique((Idents(dat.sub)))),invert=TRUE)]
ids <- ids[!is.na(ids)]
dat.subNo12 <- subset(dat.sub,idents = ids)

# get motifs 
# add motif information
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
dat.subNo12 <- AddMotifs(object = dat.subNo12,assay = "peaks",genome = BSgenome.Hsapiens.UCSC.hg38,pfm = pfm)

# RT vs no RT NTCs only
res <- FindMarkers(object = dat.subNo12,ident.1 = "non-targeting_NoRT", ident.2 = "non-targeting_5x1.8",min.pct = 0.05,test.use = 'LR')
res$gene <- annot[nearest(StringToGRanges(row.names(res)), annot, ignore.strand=TRUE)]$gene_name
res <- res[order(res$avg_log2FC),]
top.da.peak <- rownames(res[res$p_val < 0.05, ])

open.peaks <- AccessiblePeaks(dat.subNo12, idents = c("non-targeting_NoRT", "non-targeting_5x1.8"))

# match the overall GC content in the peak set
meta.feature <- GetAssayData(dat.subNo12, assay = "peaks", slot = "meta.features")
peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[top.da.peak, ],
  n = 5000
)

enriched.motifs <- FindMotifs(
  object = dat.subNo12,assay = "peaks",
  features = top.da.peak
)

MotifPlot(
  object = dat.subNo12,assay = "peaks",
  motifs = head(rownames(enriched.motifs))
)

dat.subNo12 <- RunChromVAR(
  object = dat.subNo12,assay = "peaks",
  genome = BSgenome.Hsapiens.UCSC.hg38
)

DefaultAssay(dat.subNo12) <- 'chromvar'

# look at the activity of EGR1
p1 <- DimPlot(dat.subNo12, label = TRUE, pt.size = 0.1)

p2 <- FeaturePlot(
  object = dat.subNo12,
  features = "MA0162.4",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p1 + p2

# get violin plot data for EGR1 
VlnPlot(dat.subNo12,features = "MA0162.4") + NoLegend()

# tsne of chromvar data 
dat.subNo12 <- RunTSNE(dat.subNo12, assay = 'chromvar',seed.use = 1,tsne.method = "Rtsne",dim.embed = 2)
DimPlot(dat.subNo12, label = TRUE, repel = TRUE, reduction = "tsne")

# get summary chromvar
MeanChromVars <- AverageExpression(dat.subNo12,assays = "chromvar")
MeanChromVars <- MeanChromVars$chromvar
row.names(MeanChromVars) <- enriched.motifs[row.names(MeanChromVars),"motif.name"]

write.table(MeanChromVars, 'MeanChromVars.csv', sep='\t',quote=FALSE, col.names = NA)

# prepare chromvar heatmap with metadata 
expressedGenes <- as.character(read.table('/raleighlab/data1/liuj/schwannoma/scrna_cells/analysis/expressedGenes_in_HEI193_max10.txt',header=TRUE)[,"x"])
screenTbl <- read.table('../epigHitsTbl_enrichments_forlibrary_withGO.csv', sep=',',row.names=1,header=TRUE)

metaData <- data.frame(row.names = colnames(MeanChromVars), id = colnames(MeanChromVars))   
metaData <- data.frame(metaData, gene = str_split_fixed(metaData$id, "_", 2)[,1], cond = str_split_fixed(metaData$id, "_", 2)[,2])
metaData <- data.frame(metaData, GO = screenTbl[as.character(metaData$gene),"GO"], gamma = screenTbl[as.character(metaData$gene),"gamma_avg"], rho = screenTbl[as.character(metaData$gene),"rho_avg"])
metaData[metaData$gene == "EP300","GO"] <- "Histone acetyltransferase"
metaData[metaData$gene == "BRD4","GO"] <- "BRD4"
metaData[metaData$gene == "CTCF","GO"] <- "CTCF"
metaData[metaData$gene == "SMARCB1","GO"] <- "SWI SNF"
metaData[metaData$gene == "non-targeting","GO"] <- "non-targeting"
metaData$GOcond <- paste(metaData$GO, metaData$cond,sep="_")

metaDataMotifs <- enriched.motifs
row.names(metaDataMotifs) <- metaDataMotifs$motif.name
metaDataMotifs$nlog2pval <- -log(metaDataMotifs$p.adjust,2)

# subset the data and rank by enriched motifs 
MeanChromVarsSub <- MeanChromVars[,grep("12$",invert=TRUE,colnames(MeanChromVars))]
MeanChromVarsSub <- MeanChromVarsSub[intersect(row.names(MeanChromVarsSub),expressedGenes),]


scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

zscores <- scale_rows(MeanChromVarsSub)

# plot it 
png("heatmap_motif_chromvar_expressed.png",width = 16000, height=16000,res = 600)
pheatmap(log(MeanChromVarsSub+1,2),border_color = NA,cellwidth = 10,cellheight = 10,breaks = seq(-5,5,10/50),annotation_col = metaData[,c("GO","cond","gamma","rho")], annotation_row = data.frame(row.names=row.names(metaDataMotifs),metaDataMotifs[,c("nlog2pval")]),
         color=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(50))
dev.off()

png("heatmap_motif_chromvar_expressed_zscore.png",width = 16000, height=16000,res = 600)
pheatmap(zscores,border_color = NA,cellwidth = 10,cellheight = 10,breaks = seq(-2,2,4/50),annotation_col = metaData[,c("GO","cond","gamma","rho")], annotation_row = data.frame(row.names=row.names(metaDataMotifs),metaDataMotifs[,c("nlog2pval")]),
         color=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(50))
dev.off()

# export table for exploration 

MeanChromVarsSubMeta <- rbind(t(metaData[colnames(MeanChromVarsSub),c("GO","cond","gamma","rho")]),MeanChromVarsSub)
MeanChromVarsSubMeta <- data.frame(nlog2pval = metaDataMotifs[row.names(MeanChromVarsSubMeta),"nlog2pval"], motif.enrich = metaDataMotifs[row.names(MeanChromVarsSubMeta),"fold.enrichment"], MeanChromVarsSubMeta)
write.table(MeanChromVarsSubMeta, 'MeanChromVars_metadata.csv', sep='\t',quote=FALSE, col.names = NA)

# further analysis of SAGA complex factors at EGR1/TP53
dat.subNo12$GOcond <- dat.subNo12$sgRNACond
dat.subNo12$GOcond <- metaData[as.character(dat.subNo12$GOcond),"GOcond"]

Idents(dat.subNo12) <- dat.subNo12$GOcond
DefaultAssay(dat.subNo12) <- 'chromvar'
VlnPlot(dat.subNo12,features = "MA0106.3") + NoLegend()
ggsave(filename = "vs_snarc_chromvar_violin_TP53_motifs.pdf",width = 8, height=7)

########################################################
# get p53 data
x <- scan("/scratch/liuj/genomes/ENCODE_TF_ChIP-seq_2015.txt", what="", sep="\n")
y <- strsplit(x, "[\t]")
names(y) <- sapply(y, `[[`, 1)
y <- lapply(y, `[`, -1)
y <- lapply(y, `[`, -1)
encodeChIP <- y

ARCHS4 <- read.table('/scratch/liuj/genomes/ARCHS4_TFs_Coexp.txt',sep='\t')
ARCHS4 <- ARCHS4[,-2]
ARCHS4$V1 <- gsub(" human tf ARCHS4 coexpression","",ARCHS4$V1)
row.names(ARCHS4) <- as.character(ARCHS4$V1)
ARCHS4 <- ARCHS4[,-1]

# get DARS for gene of interest and compared to NTC RT #
smp <- "TADA1_5x1.8"
smp2 <- "TADA1_NoRT"
res <- data.frame(resListAll[smp])
res$gene <- annot[nearest(StringToGRanges(row.names(res)), annot, ignore.strand=TRUE)]$gene_name
res <- res[order(res[,1]),]
res <- res[order(res[,2]),]

resSub <- res[grep(paste(as.character(ARCHS4["TP53",]),collapse="$|^"),res$gene),] 
resSub <- resSub[resSub[,2] < 0,]

# get all KD enriched DARS
resSub <- res[res[,2] < 0,]

# try rGreat - use background as all open peaks in NTC no RT + RT conditions + DARS in condition of interest - look at enriched in KD only
job = submitGreatJob(StringToGRanges(row.names(res[res[,1] < 0.05 & res[,2] > 0,])), species = "hg38", bg = c(StringToGRanges(open.peaks),StringToGRanges(row.names(res))))
tb = getEnrichmentTables(job)

plotRegionGeneAssociationGraphs(job)

# extract GREAT terms 
n = 2 # which GO terms? 
dfPlot <- tb[[n]][tb[[n]]$Hyper_Raw_PValue < 0.05 & tb[[n]]$Hyper_Fold_Enrichment > 5,]
dfPlot <- dfPlot[,c("ID","name","Hyper_Fold_Enrichment","Hyper_Raw_PValue")]
dfPlot

df <- plotRegionGeneAssociationGraphs(job, ontology = "GO Biological Process", termID = "GO:0021535") # GO term of choice to extract 
df

# look at linkages
goi <- "TP53"
Idents(dat.subNo12) <- "sgRNACond"
dat.subPlot <- subset(dat.subNo12,idents = c("non-targeting_NoRT","non-targeting_5x1.8",smp,smp2))
dat.subPlot <- LinkPeaks(object = dat.subPlot, peak.assay = "peaks", expression.assay = "SCT", genes.use = goi)

DefaultAssay(dat.subPlot) <- "ATAC"
CoveragePlot(dat.subPlot, region = goi,features = goi,extend.upstream = 10000,extend.downstream = 10000)

# look at average profile plots at motifs
motifsdf <- GetMotifData(dat.subNo12,assay = 'peaks')
colnames(motifsdf) <- ConvertMotifID(object = dat.subNo12, id = colnames(motifsdf), assay='peaks') 

# get overlap motif regions with DARS
motifoi <- motifsdf[,"TP53"]
motifoi <- motifoi[motifoi > 0]

x <- findOverlaps(StringToGRanges(names(motifoi)),StringToGRanges(row.names(resSub)),minoverlap = 1) #overlap between motif of interest and also DAR for given gene
x <- motifoi[queryHits(x)]
x <- StringToGRanges(names(x))
x <- resize(x = x,width = 100)

dat.subNo12 <- Footprint(object = dat.subNo12, genome = BSgenome.Hsapiens.UCSC.hg38,assay = "ATAC",regions = x, key = "TP53_DARS", in.peaks =TRUE)

p1 <- PlotFootprint(dat.subNo12, features = "TP53_DARS",assay = "peaks",idents = c("non-targeting_NoRT","non-targeting_5x1.8",smp,smp2))
p1 + patchwork::plot_layout(ncol = 1)

ggsave("snARC_footprint_TP53.png", width = 6, height = 6)

DefaultAssay(dat.subNo12) <- "chromvar"
FeaturePlot(object = dat.subNo12, features = "MA1513.1",min.cutoff = 'q10', max.cutoff = 'q90', pt.size = 0.1)

# cicero
# peak - RNA linkage analysis
# coverage plot of average genes - use deeptools of subsetted bam files or Granges 
# distribution of distances to TSS for DARs
# differential TSS enrichment?

# generate matrix of all DARs

resListAll <- resListAll.LR

# get all DE genes 
DARSAll <- unique(unlist(lapply(resListAll,function(x) row.names(x[as.numeric(x$p_val) < 0.05 & abs(as.numeric(x$avg_log2FC)) > 0.25,]))))
DARSAll <- DARSAll[!is.na(DARSAll)]
DARSAll <- DARSAll[grep("NA",DARSAll,invert=TRUE)]

# get distances to genes 
DARSAll.Dist <- data.frame(region = DARSAll, dist = rep(0,length(DARSAll)), gene = rep(NA, length(DARSAll)))
DARSAll.Dist$gene <-  annot[nearest(StringToGRanges(DARSAll.Dist$region), annot, ignore.strand=TRUE)]$gene_name
DARSAll.Dist$dist <- distanceToNearest(StringToGRanges(DARSAll.Dist$region), annot, ignore.strand=TRUE)@elementMetadata$distance
row.names(DARSAll.Dist) <- DARSAll.Dist$region

# build log2FC matrix 
#resListAll <- resListAll.MAST
#resListAll <- resListAll[names(resListAll) != "ENY2_NoRT"] # empty for the MAST output

for (i in 1:length(resListAll)){
  resListAll[[i]] <- data.frame(gene = row.names(resListAll[[i]]),resListAll[[i]])
}

deMat <- ldply(lapply(resListAll,function(x) x[,c("gene","avg_log2FC")]),data.frame)
deMat <- reshape(deMat, idvar = "gene", v.names = "avg_log2FC",timevar = ".id", direction = "wide")
deMat <- deMat[!is.na(deMat$gene),]
row.names(deMat) <- deMat$gene
colnames(deMat) <- gsub("avg_log2FC.","",colnames(deMat))
deMat <- deMat[,-1]
deMat <- deMat[intersect(DARSAll,row.names(deMat)),]
deMat[is.na(deMat)] <- 0
deMat[deMat == Inf] <- 0
deMat[deMat == -Inf] <- 0

  png("snARC_DARs_LRtest_sig.png",width = 12000, height=6000,res = 600)
  pheatmap(-deMat,scale = "row",border_color = NA,annotation_col = metaData[,c("GO","cond","gamma","rho")], 
           #annotation_row = data.frame(row.names=row.names(DARSAll.Dist), log2dist = log(DARSAll.Dist[,c("dist")]+1,2)),
           annotation_names_row = FALSE, breaks = seq(-1,1,(2/50)), color=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(50))
  dev.off()

# export annotated heatmap for exploration 
deMatMeta <- rbind(t(metaData[colnames(deMat),c("GO","cond","gamma","rho")]),-deMat)
deMatMeta <- data.frame(closest_gene = DARSAll.Dist[row.names(deMatMeta),"gene"], log2dist = log(DARSAll.Dist[row.names(deMatMeta),"dist"]+1,2), deMatMeta)
write.table(deMatMeta, 'DARMat_MAST_with_metadata_sig.csv', sep='\t',quote=FALSE, col.names = NA)

# Boxplot of avg_log2fc of genes by gene and GO term 
# Violin plot of of distances between DARs and genes by GO term and by gene
# get better p53 and EGR1 peaks for individual gene analysis 
dfLog2fc <- ldply(resListAll, function(x) x[as.numeric(x$p_val) < 0.05 & abs(as.numeric(x$avg_log2FC)) > 0.25,])
dfLog2fc$gene <- str_split_fixed(dfLog2fc$.id, "_", 2)[,1]
dfLog2fc$cond <- str_split_fixed(dfLog2fc$.id, "_", 2)[,2]
dfLog2fc$GO <- metaData[as.character(dfLog2fc$.id),"GO"]
dfLog2fc$GOcond <- metaData[as.character(dfLog2fc$.id),"GOcond"]
dfLog2fc$dist <- distanceToNearest(StringToGRanges(dfLog2fc$gene.1), annot, ignore.strand=TRUE)@elementMetadata$distance
dfLog2fc$closestGene <- annot[nearest(StringToGRanges(dfLog2fc$gene.1), annot, ignore.strand=TRUE)]$gene_name
dfLog2fc$distCategory <-  dfLog2fc$dist
dfLog2fc[dfLog2fc$dist < 1000,"distCategory"] <- "Promoter"
dfLog2fc[dfLog2fc$dist >= 1000 & dfLog2fc$dist < 10000,"distCategory"] <- "Enhancer < 10 kb"
dfLog2fc[dfLog2fc$dist >= 10000,"distCategory"] <- "Enhancer >= 10 kb"

plot <- ggplot(dfLog2fc) + 
  geom_boxplot(aes(x=gene, y=-avg_log2FC, fill=cond),size=0.25,outlier.size = 0.25) +
  #scale_x_continuous(limits = c(0,10),breaks = seq(0,10,1),expand=c(0.01,0.01)) +
  #scale_y_continuous(expand=c(0.01,0.01), limits = c(0,11000)) +
  labs(title = "DARS normalized to non-targeting within RT conditions") +
  #facet_wrap(~.id, scales="free", ncol = 2) +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_blank(),panel.background=element_blank(),axis.line = element_line(size=0.25))
ggsave("snarc_DARs_boxplot_genes.png", plot, width =12, height = 5)

plot <- ggplot(dfLog2fc) + 
  geom_boxplot(aes(x=GO, y=-avg_log2FC, fill=cond),size=0.25,outlier.size = 0.25) +
  #scale_x_continuous(limits = c(0,10),breaks = seq(0,10,1),expand=c(0.01,0.01)) +
  #scale_y_continuous(expand=c(0.01,0.01), limits = c(0,11000)) +
  labs(title = "DARS normalized to non-targeting within RT conditions") +
  #facet_wrap(~.id, scales="free", ncol = 2) +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_blank(),panel.background=element_blank(),axis.line = element_line(size=0.25))
ggsave("snarc_DARs_boxplot_GO.png", plot, width =7, height = 4)

plot <- ggplot(dfLog2fc) + 
  geom_boxplot(aes(x=gene, y=log(dist+1,2), fill=cond),size=0.25,outlier.size = 0.25) +
  #scale_x_continuous(limits = c(0,10),breaks = seq(0,10,1),expand=c(0.01,0.01)) +
  #scale_y_continuous(expand=c(0.01,0.01), limits = c(0,11000)) +
  labs(title = "Distance between DARS and closest gene") +
  #facet_wrap(~.id, scales="free", ncol = 2) +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_blank(),panel.background=element_blank(),axis.line = element_line(size=0.25))
ggsave("snarc_DARs_boxplot_distanceToNearestGene.png", plot,  width =12, height = 5)

# distance plots - category
distComp <- data.frame(table(dfLog2fc[,c(".id","distCategory")]))
distComp$gene <- str_split_fixed(distComp$.id, "_", 2)[,1]
distComp$cond <- str_split_fixed(distComp$.id, "_", 2)[,2]
distComp$GO <- metaData[as.character(distComp$.id),"GO"]
distComp$GOcond <- metaData[as.character(distComp$.id),"GOcond"]

plot <- ggplot(distComp) +
  geom_bar(aes(x = .id, y = Freq, fill = distCategory), position = "fill",stat = 'identity') +
  scale_y_continuous(labels = scales::percent,expand = c(0.005,0.005)) +
  labs(title = "Distance between DARs and nearest gene", x="", y="Percentage") +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=45,hjust=1,vjust=1),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_rect(colour='black'))
ggsave("snarc_DARS_distance_categorical.png", plot, width = 12, height = 4)

plot <- ggplot(distComp) +
  geom_bar(aes(x = GOcond, y = Freq, fill = distCategory), position = "fill",stat = 'identity') +
  scale_y_continuous(labels = scales::percent,expand = c(0.005,0.005)) +
  labs(title = "Distance between DARs and nearest gene", x="", y="Percentage") +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=45,hjust=1,vjust=1),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_rect(colour='black'))
ggsave("snarc_DARS_distance_categorical_GO.png", plot, width = 7, height = 4)

# Just hone in on the important motifs (top 100 enriched) whose TFs are also expressed 
MeanChromVarsSubEnriched <- MeanChromVarsSub[intersect(row.names(MeanChromVarsSub),toupper(enriched.motifs[1:100,"motif.name"])),]

# plot it 
png("heatmap_motif_chromvar_expressed_enriched100.png",width = 14000, height=6000,res = 600)
pheatmap(log(MeanChromVarsSubEnriched+1,2),border_color = NA,cellwidth = 10,cellheight = 10,breaks = seq(-1,1,2/50),annotation_col = metaData[,c("GO","cond","gamma","rho")], annotation_row = data.frame(row.names=row.names(metaDataMotifs),metaDataMotifs[,c("nlog2pval")]),
         color=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(50))
dev.off()

# calculate different between noRT and RT conditions 
MeanChromVarsSubEnrichedDiff <- data.frame(row.names = row.names(MeanChromVarsSubEnriched), matrix(data = 0,nrow = nrow(MeanChromVarsSubEnriched),ncol = ncol(MeanChromVarsSubEnriched)/2))
colnames(MeanChromVarsSubEnrichedDiff) <- unique(str_split_fixed(colnames(MeanChromVarsSubEnriched), "_", 2)[,1])

for (i in 1:nrow(MeanChromVarsSubEnrichedDiff)){
  for (j in 1:ncol(MeanChromVarsSubEnrichedDiff)){
    motif = row.names(MeanChromVarsSubEnrichedDiff)[i]
    perturbation = colnames(MeanChromVarsSubEnrichedDiff)[j]
    paste(perturbation,"NoRT",sep="_")
    paste(perturbation,"5x1.8",sep="_")
    MeanChromVarsSubEnrichedDiff[i,j] <- MeanChromVarsSubEnriched[motif,paste(perturbation,"5x1.8",sep="_")] - MeanChromVarsSubEnriched[motif,paste(perturbation,"NoRT",sep="_")] 
  }
}

MeanChromVarsSubEnrichedDiff <- MeanChromVarsSubEnrichedDiff[row.names(metaDataMotifs[order(-metaDataMotifs$nlog2pval),]),] # order by enrichment
MeanChromVarsSubEnrichedDiff <- MeanChromVarsSubEnrichedDiff[grep("NA",row.names(MeanChromVarsSubEnrichedDiff),invert=TRUE),]
MeanChromVarsSubEnrichedDiff <- MeanChromVarsSubEnrichedDiff[,colnames(MeanChromVarsSubEnrichedDiff[,order(-as.numeric(MeanChromVarsSubEnrichedDiff["EGR1",]))])]

pdf("heatmap_motif_chromvar_expressed_enriched100_RTdiff_clustered.pdf",width = 18, height=12)
pheatmap(MeanChromVarsSubEnrichedDiff,border_color = NA,cellwidth = 10,cellheight = 10,breaks = seq(-3,3,6/50), annotation_row = data.frame(row.names=row.names(metaDataMotifs),metaDataMotifs[,c("nlog2pval")]),
         cluster_rows = TRUE, cluster_cols = TRUE,color=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(50))
dev.off()

###########################################
# Reformat heatmap using complex heatmaps #

metaDataChromVar <- metaData
row.names(metaDataChromVar) <- make.unique(as.character(metaDataChromVar$gene))
metaDataChromVar["non-targeting","rho"] <- mean(0.00452918636789579, -0.00221716500559651)
metaDataChromVar["non-targeting","gamma"] <- mean(0.00724089117087796, 0.000262435466018865)
metaDataChromVar <- metaDataChromVar[!is.na(metaDataChromVar$rho),]
metaDataChromVar <- metaDataChromVar[intersect(row.names(metaDataChromVar),colnames(MeanChromVarsSubEnrichedDiff)),]

# remove additions
metaDataChromVar <- metaDataChromVar[c(screenTbl[screenTbl$thresh.rho > 4,"gene"],"non-targeting"),]
metaDataChromVar <- metaDataChromVar[grep("NF2",invert=TRUE,row.names(metaDataChromVar)),]
metaDataChromVar$GO <- factor(as.character(metaDataChromVar$GO), levels = unique(metaDataChromVar$GO))

pdf("heatmap_motif_chromvar_expressed_enriched100_RTdiff_clustered_complex_rhoonly.pdf",width = 10, height=8)
Heatmap(as.matrix(MeanChromVarsSubEnrichedDiff[,intersect(row.names(metaDataChromVar),colnames(MeanChromVarsSubEnrichedDiff))]),
        name = "Motif Deviation\n(RT - 0Gy)",
        top_annotation = HeatmapAnnotation(GO = metaDataChromVar$GO, 
                                           col = list(GO = c("Histone demethylase" = "#9E0142", "Histone acetyltransferase" = "#D53E4F",
                                                              "Histone deacetylase" = "#F46D43", "SAGA complex" = "#FDAE61",
                                                              "HBO1 complex" = "#FEE08B", "Histone methyltransferase" = "#FFFFBF",
                                                              "NuA4 complex" = "#E6F598", "NuRD complex" = "#ABDDA4","PRC2 complex" = "#66C2A5",
                                                              "SWI SNF" = "#3288BD", "Glycosylase" = "#5E4FA2", "non-targeting" = "gray50")),
                                           gamma = anno_barplot(metaDataChromVar$gamma, bar_width = 0.9,ylim = c(-1,1)), 
                                           rho = anno_barplot(metaDataChromVar$rho, bar_width = 0.9)),
        #column_split = factor(sgRNAClassRT$cond, levels = c("0Gy","9Gy","12_5Gy")),
        cluster_column_slice = FALSE,
        width = 30*unit(4, "mm"),
        height = 27*unit(4, "mm"),
        col = colorRamp2(seq(-3, 3, length = 3), c("#2268AD", "#EEEEEE", "#B31B2C"))
)
dev.off()


# Go back to DARS associated with ASH2L, TADA1/2B, KDM1A/5C at EGR1 genes. Make metagene tables
# get EGR1 target genes and subset the deMat table 
EGR1Targets <- read.table('/scratch/liuj/genomes/EGR1_01.v7.5.1.grp',sep='\t',skip=1)
EGR1Targets <- as.character(EGR1Targets$V1)

DARSAll.Dist.EGR1 <- DARSAll.Dist[grep(paste(EGR1Targets,collapse="$|^"),DARSAll.Dist$gene),]
deMatEGR1Targets <- deMat[row.names(DARSAll.Dist.EGR1),]

png("snARC_DARs_LRtest_sig_EGR1targets.png",width = 12000, height=3000,res = 600)
pheatmap(-deMatEGR1Targets,scale = "row",border_color = NA,annotation_col = metaData[,c("GO","cond","gamma","rho")], 
         #annotation_row = data.frame(row.names=row.names(DARSAll.Dist), log2dist = log(DARSAll.Dist[,c("dist")]+1,2)),
         annotation_names_row = FALSE, breaks = seq(-1,1,(2/50)), color=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(50))
dev.off()

df <- deMatEGR1Targets
df[df == 0] <- NA

# generate boxplot summaries of DARs at EGR1 sites 
# average profile plots of ATAC data 
smp <- "KDM5C_NoRT"
smp2 <-"KDM5C_5x1.8"
dars <- resListAll[[smp2]]
dars <- dars[dars$p_val < 0.05,]
dars$gene.1 <- DARSAll.Dist[row.names(dars),"gene"]
dars <- dars[grep(paste(EGR1Targets,collapse="$|^"),dars$gene.1),]

regionsOI <- StringToGRanges(row.names(dars))
strand(regionsOI) <- "+"

df <- TSSEnrichment(dat.subNo12, fast=FALSE, assay = 'ATAC',verbose=TRUE,tss.positions = regionsOI,)
TSSPlot(df,assay = "ATAC",idents = c(smp,smp2,"non-targeting_NoRT","non-targeting_5x1.8"))

goi <- "NAB2"
Idents(dat.subNo12) <- "sgRNACond"
dat.subPlot <- subset(dat.subNo12,idents = c("non-targeting_NoRT","non-targeting_5x1.8",smp, smp2))
dat.subPlot <- LinkPeaks(object = dat.subPlot, peak.assay = "peaks", expression.assay = "SCT", genes.use = goi)

DefaultAssay(dat.subPlot) <- "ATAC"
CoveragePlot(dat.subPlot, region = goi,features = goi,extend.upstream = 10000,extend.downstream = 10000)

# Export RDS so Jason can run Cicero
save(dat.subNo12,file="dat.subNo12.rds")

## Analysis of other DARs from motif heatmap 
x <- scan("/scratch/liuj/genomes/ENCODE_TF_ChIP-seq_2015.txt", what="", sep="\n")
y <- strsplit(x, "[\t]")
names(y) <- sapply(y, `[[`, 1)
y <- lapply(y, `[`, -1)
y <- lapply(y, `[`, -1)
encodeChIP <- y

# load custom gene sets
gset <- read.table('/scratch/liuj/genomes/jaspar_KLF4.txt',sep='\t')
gset <- as.character(gset$V1)

# look at TCF3 genes in the setting of TADA1 KD 
smp <-"SETDB1_NoRT"
smp2 <-"SETDB1_5x1.8"
dars <- resListAll[[smp2]]
dars <- dars[dars$p_val < 0.05,]
dars$gene.1 <- DARSAll.Dist[row.names(dars),"gene"]
dars <- dars[grep(paste(encodeChIP[["TCF3 GM12878 hg19"]],collapse="$|^"),dars$gene.1),]
#dars <- dars[grep(paste(as.character(ARCHS4["HES1",]),collapse="$|^"),dars$gene.1),]

regionsOI <- StringToGRanges(row.names(dars))
strand(regionsOI) <- "+"

df <- TSSEnrichment(dat.subNo12, fast=FALSE, assay = 'ATAC',verbose=TRUE,tss.positions = regionsOI)
TSSPlot(df,assay = "ATAC",idents = c(smp,smp2,"non-targeting_NoRT","non-targeting_5x1.8"))

goi <- "TGFBR1"
Idents(dat.subNo12) <- "sgRNACond"
dat.subPlot <- subset(dat.subNo12,idents = c("non-targeting_NoRT","non-targeting_5x1.8",smp, smp2))
dat.subPlot <- LinkPeaks(object = dat.subPlot, peak.assay = "peaks", expression.assay = "SCT", genes.use = goi)

DefaultAssay(dat.subPlot) <- "ATAC"
CoveragePlot(dat.subPlot, region = goi,features = goi,extend.upstream = 10000,extend.downstream = 10000)

# 8/12/22
# load Jason's cicero output
dat.subNo12Cicero <- readRDS('../dat.subNo12_chromatin_contacts.rds')

# load KLF Encode ChIP data
gset <- read.table('../KLF6_encode_chip.csv',sep=',',header=TRUE) # need to fix 
gset=apply(gset[,1:3],1,function(x) paste(x,collapse = "-",sep=''))
gset=gsub(" ", "", gset, fixed = TRUE)
gset=gset[grep("chr1|chr2|chr3|chr4|chr5|chr6|chr7|chr8|chr9|chr10|chr11|chr12|chr13|chr14|chr15|chr16|chr17|chr18|chr19|chr20|chr21|chr22|chrX|chrY",gset)]
gset=gset[grep("random",gset,invert = TRUE)]
gset <- unique(annot[nearest(StringToGRanges(gset), annot, ignore.strand=TRUE)]$gene_name)

# look at ChIP target genes in the setting of KDM1A/5C KD 
smp <-"KDM5C_NoRT"
smp2 <-"KDM5C_5x1.8"
smp3 <- "KDM1A_NoRT"
smp4 <- "KDM1A_5x1.8"
dars <- resListAll[[smp2]]
#dars <- dars[dars$p_val < 0.1,] # all detected sites in a given perturbation
dars$gene.1 <- DARSAll.Dist[row.names(dars),"gene"]
dars <- dars[match(gset,dars$gene.1),]
dars <- dars[!is.na(dars$gene),]

regionsOI <- StringToGRanges(row.names(dars))
strand(regionsOI) <- "+"

df <- TSSEnrichment(dat.subNo12Cicero, fast=FALSE, assay = 'ATAC',verbose=TRUE,tss.positions = regionsOI)
TSSPlot(df,assay = "ATAC",idents = c(smp,smp2,"non-targeting_NoRT","non-targeting_5x1.8"))

goi <- "TGFBR1"
Idents(dat.subNo12) <- "sgRNACond"
dat.subPlot <- subset(dat.subNo12,idents = c("non-targeting_NoRT","non-targeting_5x1.8",smp, smp2))
dat.subPlot <- LinkPeaks(object = dat.subPlot, peak.assay = "peaks", expression.assay = "SCT", genes.use = goi)

DefaultAssay(dat.subPlot) <- "ATAC"
CoveragePlot(dat.subPlot, region = goi,features = goi,extend.upstream = 10000,extend.downstream = 10000)

##########################################################
# Standardize for any ChIP-seq data and any perturbation #
# analysis of ChIP-seq from ENCODE 
#gset <- read.table('../KLF13_chip_bedidr_threshold_peakENCFF453MMH.bed',sep='\t',header=FALSE)
gset <- read.table('../TCF3_IDR_threshold_peaks.csv',sep=',',header=TRUE)

gset=apply(gset[,1:3],1,function(x) paste(x,collapse = "-",sep=''))
gset <- gsub(' ','',gset)
gset=gset[grep("chr1|chr2|chr3|chr4|chr5|chr6|chr7|chr8|chr9|chr10|chr11|chr12|chr13|chr14|chr15|chr16|chr17|chr18|chr19|chr20|chr21|chr22|chrX|chrY",gset)]
gset=gset[grep("random",gset,invert = TRUE)]
gset <- unique(annot[nearest(StringToGRanges(gset), annot, ignore.strand=TRUE)]$gene_name)

# process plots with motif and perturbation of interest
motif <- "TCF3"
goi <- "SETDB1"

smp <- paste(goi,"_NoRT",sep='')
smp2 <- paste(goi,"_5x1.8",sep='')
dars <- resListAll[[smp2]]
#dars <- dars[dars$p_val < 0.1,] # all detected sites in a given perturbation
dars$gene.1 <- DARSAll.Dist[row.names(dars),"gene"]
dars <- dars[match(gset,dars$gene.1),] # Match ChIP-seq data with DARS
dars <- dars[!is.na(dars$gene),]
regionsOI <- StringToGRanges(row.names(dars))
regionsOI <- resize(regionsOI, width = 1, fix = "center") # center around the DAR
strand(regionsOI) <- "+"

# Export DARS for others 
#write.table(dars, paste(goi,"_perturb_at_",motif,".csv",sep=''), sep='\t',quote=FALSE, row.names=FALSE)

# calcuate enrichment values
df <- TSSEnrichment(dat.subNo12Cicero, fast=FALSE, assay = 'ATAC',verbose=TRUE,tss.positions = regionsOI)
plot <- TSSPlot(df,assay = "ATAC",idents = c(smp,smp2,"non-targeting_NoRT","non-targeting_5x1.8"))

# try to export TSSenrichment data and plot them together 
plot <- plot$data
plot$group <- factor(plot$group,levels = c(smp2,smp,"non-targeting_5x1.8","non-targeting_NoRT"))

ggplot(plot) + 
  #geom_line(aes(x=position, y=norm.value, color=group),size=0.25) +
  geom_smooth(method='loess',aes(x=position, y=norm.value, color=group),size=1,span=0.5) +
  labs(title = paste(goi,"perturbation at",motif,"motifs")) +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_blank(),panel.background=element_blank(),axis.line = element_line(size=0.25))
ggsave(paste("snarc_",goi,"_pertub_",motif,"_motif_avgprofile.pdf",sep=''), width =4, height = 2.5)

# Generate violin plot for scRNA-seq   
datPlot <- subset(dat.subNo12Cicero,idents = c(smp,smp2,"non-targeting_NoRT","non-targeting_5x1.8"))
DefaultAssay(datPlot) <- "RNA"
datPlot <- NormalizeData(datPlot, normalization.method = "LogNormalize", scale.factor = 10000)
datPlot <- FindVariableFeatures(datPlot, selection.method = "vst", nfeatures = 3000)
df <- GetAssayData(datPlot,slot = 'data',assay = 'RNA')
df <- df[intersect(row.names(df),gset),]
df <- df[apply(df,1,function(x) sum(x > 0)) > 5,] # genes with at least 5 cells 

df <- data.frame(cells = colnames(df), avg = apply(df,2,function(x) mean(x)))
df$sgRNACond <- as.character(Idents(dat.subNo12Cicero)[row.names(df)])

ggplot(df) + 
  geom_boxplot(aes(x=sgRNACond, y=avg, fill=sgRNACond),size=0.25,outlier.size = 1) +
  #scale_x_continuous(limits = c(0,10),breaks = seq(0,10,1),expand=c(0.01,0.01)) +
  #scale_y_continuous(expand=c(0.01,0.01), limits = c(0,11000)) +
  labs(title = paste("snRNA",goi,"at",motif,"motifs"),y = "Average Expression") +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), legend.position = "none",
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_blank(),panel.background=element_blank(),axis.line = element_line(size=0.25))
ggsave(paste("snRNA_",goi,"_pertub_",motif,"_motif_boxplot.pdf",sep=''), width = 2, height = 3)

write.table(df, paste("snRNA_",goi,"_pertub_",motif,"_motif_boxplot.txt",sep=''), sep='\t',quote=FALSE, col.names = NA)

# Try addmodule score for RNA 
datPlot <-AddModuleScore(datPlot,features = list(intersect(VariableFeatures(datPlot),gset)), name = 'module')
VlnPlot(datPlot, 'module')

# plot sample genes 
goi <- "STX2"
Idents(dat.subNo12Cicero) <- "sgRNACond"
dat.subPlot <- subset(dat.subNo12Cicero,idents = c("non-targeting_NoRT","non-targeting_5x1.8",smp, smp2))
dat.subPlot <- LinkPeaks(object = dat.subPlot, peak.assay = "peaks", expression.assay = "SCT", genes.use = goi)

DefaultAssay(dat.subPlot) <- "ATAC"


# export table for exploration 

MeanChromVarsSubMeta <- rbind(t(metaData[colnames(MeanChromVarsSub),c("GO","cond","gamma","rho")]),MeanChromVarsSub)
MeanChromVarsSubMeta <- data.frame(nlog2pval = metaDataMotifs[row.names(MeanChromVarsSubMeta),"nlog2pval"], motif.enrich = metaDataMotifs[row.names(MeanChromVarsSubMeta),"fold.enrichment"], MeanChromVarsSubMeta)
write.table(MeanChromVarsSubMeta, 'MeanChromVars_metadata.csv', sep='\t',quote=FALSE, col.names = NA)
