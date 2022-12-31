## Analysis of clinical single cell RNA-seq with integration of single nuclei RNA-seq of VS
## Only keep samples with > 100 cells 
library(Seurat)
library(ggplot2)
library(ggthemes)
library(GGally)
library(dplyr)
library(harmony)
library(reshape)
library(data.table)
library(WriteXLS)
library(sctransform)
library(RColorBrewer)
library(clusterProfiler)
library(msigdbr)
library(plyr)
library(vegan)
library(SCpubr)

setwd('/raleighlab/data1/liuj/schwannoma/scrna/analysis_scrna')

# read in data
smpid <- dir(file.path("/raleighlab/data1/liuj/schwannoma/scrna/cellranger_out/"))
smpidShort <- smpid
smpidShort <- gsub("sn_VS_","",smpidShort)
smpidShort <- gsub("vs_sn_","",smpidShort)
smpidShort <- gsub("VS_","",smpidShort)

data.list <- list()

for (i in 1:length(smpid)){
  df <- Read10X(data.dir = paste("/media/data2/liuj/vestibular/scrna/cellranger_out/",smpid[i],"/outs/filtered_feature_bc_matrix/",sep=''))
  data.list[[i]] <- CreateSeuratObject(counts = df, project = smpid[i], min.cells = 3, min.features = 100)
}

names(data.list) <- smpidShort

dat.all <- merge(x = data.list[[1]], y = data.list[2:length(smpidShort)])

# shorten identity names
dat.all$orig.ident <- gsub("sn_VS_","",dat.all$orig.ident)
dat.all$orig.ident <- gsub("vs_sn_","",dat.all$orig.ident)
dat.all$orig.ident <- gsub("VS_","",dat.all$orig.ident)

# distinguish nuclei vs cell
dat.all$source <- dat.all$orig.ident
dat.all$source[grep("S10|S14|S15|S20|S29|S38|S6_1|S6_2|S12|S13|S24|S36|S9_1|S11|S17",dat.all$source)] <- "nuclei"
dat.all$source[grep("VSF",dat.all$source)] <- "cell"

# fraction mitochondria genes added to default qc
mito.features <- grep(pattern = "^MT-", x = rownames(x = dat.all), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = dat.all, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = dat.all, slot = 'counts'))
dat.all[['percent.mito']] <- percent.mito

# plot it 
plot <-VlnPlot(object = dat.all, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"),pt.size = 0.1,group.by = 'orig.ident')
ggsave("nc_combined_2/vsf_qc.pdf", plot, width = 12, height = 4, useDingbats = FALSE, limitsize=FALSE)
ggsave("nc_combined_2/vsf_qc.png", plot, width = 12, height = 4)

plot <- FeatureScatter(object = dat.all, feature1 = "nCount_RNA", feature2 = "percent.mito",pt.size = 0.1,group.by = 'orig.ident')
ggsave("nc_combined_2/vsf_nGenes_vs_percentmito.pdf", plot, width = 6, height = 6, useDingbats = FALSE, limitsize=FALSE)

# subset the nuclei separately with its own criteria and re-merge them 
Idents(dat.all) <- dat.all$source
dat.all.cells <- subset(x = dat.all, idents = c("cell"), subset = nFeature_RNA > 200 & percent.mito < 0.15)
dat.all.nuclei <- subset(x = dat.all, idents = c("nuclei"), subset = nFeature_RNA > 400 & percent.mito < 0.05)

dat.all.filtered <- merge(x = dat.all.cells, y = list(dat.all.nuclei))
Idents(dat.all.filtered) <- dat.all.filtered$orig.ident
dat.all.filtered <- subset(dat.all.filtered, idents = c("S13", "S24", "S36"), invert = TRUE)

# generate metadata - covariate of batch
dat.all.filtered$batch <- dat.all.filtered$orig.ident
dat.all.filtered$batch <- gsub("^S11$","r1",dat.all.filtered$batch)
dat.all.filtered$batch <- gsub("^S12$","r3",dat.all.filtered$batch)
dat.all.filtered$batch <- gsub("^S17$","r1",dat.all.filtered$batch)
dat.all.filtered$batch <- gsub("^S29$","r2",dat.all.filtered$batch)
dat.all.filtered$batch <- gsub("^S38$","r2",dat.all.filtered$batch)
dat.all.filtered$batch <- gsub("^S9_1$","r3",dat.all.filtered$batch)

# qc filters - use single cell RNA-seq filter for this 
dat.all.filtered <- SCTransform(object = dat.all.filtered, vars.to.regress = c("nCount_RNA", "percent.mito","source"), verbose = FALSE)

# linear dimensionality reduction
dat.all.filtered <- RunPCA(object = dat.all.filtered, verbose = FALSE)
plot <- DimPlot(object = dat.all.filtered, reduction = "pca",pt.size = 0.25,group.by = 'orig.ident')
ggsave("nc_combined_2/vs_nc_combined_2_pca.pdf", plot, width = 5, height = 4, useDingbats = FALSE, limitsize=FALSE)

ElbowPlot(object = dat.all.filtered,ndims = 30) 

# run harmony through seurat
dat.all.filtered <- RunHarmony(dat.all.filtered, c("source","batch"),max.iter.harmony = 20, plot_convergence = TRUE)

# Run UMAP on harmony embeddings
mindist=0.3
dat.all.filtered <- RunUMAP(dat.all.filtered, reduction = "harmony",dims = 1:20, min.dist = mindist)
plot <- DimPlot(object = dat.all.filtered, reduction = "umap", group.by = "orig.ident", pt.size=0.3, do.label=T)
ggsave(paste("nc_combined_2/vs_nc_combined_2_harmony_umap_mindist",mindist,".png",collapse="",sep=""), plot, width =8, height = 6)

plot <- DimPlot(object = dat.all.filtered, reduction = "umap", group.by = "source", pt.size=0.3, do.label=T)
ggsave(paste("nc_combined_2/vs_nc_combined_2_harmony_umap_source_mindist",mindist,".png",collapse="",sep=""), plot, width =8, height = 6)

# find clusters using harmony embeddings
dat.all.filtered <- FindNeighbors(object = dat.all.filtered, dims = 1:20, reduction = "harmony")

res=0.3
dat.all.filtered <- FindClusters(object = dat.all.filtered, reduction = "harmony", resolution = res)

# plot clusters in UMAP space 
plot <- DimPlot(object = dat.all.filtered, reduction = "umap", pt.size=0.3, do.label=T, group.by = paste("SCT_snn_res.",res,collapse="",sep=""))
ggsave(paste("nc_combined_2/vs_nc_combined_2_harmony_umap_clusters_mindist",mindist,"_res",res,".png",collapse="",sep=""), plot, width =8, height = 6)

# identify cluster markers
vsf.markers <- FindAllMarkers(object = dat.all.filtered, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25,)
vsf.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)

# marker heatmap
top10 <- vsf.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
write.table(top10,paste("nc_combined_2/top10markers__harmony_res",res,".txt",sep=''),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)

# top 50
top50 <- vsf.markers %>% group_by(cluster) %>% top_n(50, avg_logFC)
write.table(top50,paste("nc_combined_2/top50markers_harmony_res",res,".txt",sep=''),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)

# revised heatmap with normalized data
plot <- DoHeatmap(object = dat.all.filtered, features = intersect(top10$gene,row.names(dat.all.filtered@assays$SCT@scale.data)), slot = 'scale.data',draw.lines = FALSE)
ggsave(paste("nc_combined_2/vs_nc_combined_2_heatmap_harmony_top10marker_scaled_res",res,".pdf",sep=''), plot, width = 9, height = 12, useDingbats = FALSE, limitsize=FALSE)
ggsave(paste("nc_combined_2/vs_nc_combined_2_heatmap_harmony_top10marker_scaled_res",res,".png",sep=''), plot, width = 9, height = 12,limitsize = FALSE)

plot <- DoHeatmap(object = dat.all.filtered, features = intersect(top50$gene,row.names(dat.all.filtered@assays$SCT@scale.data)), slot = 'scale.data',draw.lines = FALSE)
ggsave(paste("nc_combined_2/vs_nc_combined_2_heatmap_harmony_top50marker_scaled_res",res,".pdf",sep=''), plot, width = 9, height = 24, useDingbats = FALSE, limitsize=FALSE)
ggsave(paste("nc_combined_2/vs_nc_combined_2_heatmap_harmony_top50marker_scaled_res",res,".png",sep=''), plot, width = 9, height = 24,limitsize = FALSE)

# get summary stats
clusterComp02 <- data.frame(table(dat.all.filtered[[c("orig.ident","SCT_snn_res.0.3")]]))
sampleComp02 <- data.frame(table(dat.all.filtered[[c("SCT_snn_res.0.3","orig.ident")]]))

plot <- ggplot(clusterComp02) +
  geom_bar(aes(x = SCT_snn_res.0.3, y = Freq, fill = orig.ident), position = "fill",stat = 'identity') +
  scale_y_continuous(labels = scales::percent,expand = c(0.005,0.005)) +
  labs(title = "Cluster Composition", x="Cluster", y="") +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_rect(colour='black'))
ggsave(paste("nc_combined_2/clusterComp_res",res,".pdf",sep=''), plot, width = 4, height = 4, useDingbats = FALSE, limitsize=FALSE)

plot <- ggplot(sampleComp02) +
  geom_bar(aes(x = orig.ident, y = Freq, fill = SCT_snn_res.0.3), position = "fill",stat = 'identity') +
  scale_y_continuous(labels = scales::percent,expand = c(0.005,0.005)) +
  labs(title = "Sample Composition", x="Sample", y="") +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=90,hjust=1,vjust=0),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_rect(colour='black'))
ggsave(paste("nc_combined_2/sampleComp_res",res,".pdf",sep=''), plot, width = 4, height = 4, useDingbats = FALSE, limitsize=FALSE)

# long to wide for cluster compositions
sampleCompWide <- dcast(sampleComp02, orig.ident ~ SCT_snn_res.0.3, value.var="Freq")
row.names(sampleCompWide) <- as.character(sampleCompWide$orig.ident)
sampleCompWide <- sampleCompWide[,-1]
write.table(sampleCompWide,'nc_combined_2/vs_nc_combined_2_cluster_comp.txt',sep='\t',quote=FALSE,row.names=TRUE,col.names=NA)

# Plot single cell locations on UMAP
scellEmbeddings <- read.table('vs_UMAP_res0.3_embeddings.txt',header=TRUE,row.names=1,sep='\t')
dat.all.filtered$scellLabels <- scellEmbeddings[colnames(dat.all.filtered),"cluster"]

plot <- DimPlot(object = dat.all.filtered, reduction = "umap", pt.size=0.3, do.label=T, group.by = "scellLabels")
ggsave(paste("nc_combined_2/vs_nc_combined_2_harmony_umap_scellLabeles_mindist",mindist,"_res",res,".png",collapse="",sep=""), plot, width =8, height = 6)

# export embeddings for further exploration
embeddings <- Embeddings(dat.all.filtered, reduction="umap")
embeddings <- data.frame(embeddings, combined_cluster=as.character((dat.all.filtered$SCT_snn_res.0.3[row.names(embeddings)])),
                         singleCell_cluster=as.character((dat.all.filtered$scellLabels[row.names(embeddings)])),
                         sample=as.character((dat.all.filtered$orig.ident[row.names(embeddings)])))
write.table(embeddings,'nc_combined_2/vs_nc_combined_2_UMAP_res0.3_embeddings.txt',sep='\t',quote=FALSE,row.names=TRUE,col.names=NA)

# feature plots 
# expression ordered featureplots
goi <- top50$gene
goi <- c("P2RY12")
for (i in 1:length(goi)){
  fp <- FeaturePlot(object = dat.all.filtered, reduction = "umap", features = goi[i], pt.size = 0.1)
  fp$data <- fp$data[order(fp$data[,3]),]
  ggsave(paste("nc_combined_2/featureplots/vs_nc_combined_2_umap_feature_",goi[i],".pdf",sep=''), fp, width = 3, height = 2.5, useDingbats = FALSE, limitsize=FALSE)
}

# S100B plot with smaller dots
goi <- c("S100B")
for (i in 1:length(goi)){
  fp <- FeaturePlot(object = dat.all.filtered, reduction = "umap", features = goi[i], pt.size = 0.25)
  fp$data <- fp$data[order(fp$data[,3]),]
  ggsave(paste("nc_combined_2/featureplots/vs_nc_combined_2_umap_feature_fine_",goi[i],".pdf",sep=''), fp, width = 5, height = 4, useDingbats = FALSE, limitsize=FALSE)
}

# cluster compositions 
sampleCompWideNCIE <- data.frame(type = c("IE","IE","NCE","IE","NCE","NCE","NCE","IE","IE"), sampleCompWide)
sampleCompWideNCIE <- aggregate(. ~ type, sampleCompWideNCIE, sum)
row.names(sampleCompWideNCIE) <- as.character(sampleCompWideNCIE$type)
sampleCompWideNCIE <- sampleCompWideNCIE[,-1]
sampleCompWideNCIE <- sampleCompWideNCIE / apply(sampleCompWideNCIE,1,sum)
IENCEratio <- sampleCompWideNCIE["IE",]/sampleCompWideNCIE["NCE",]
IENCEratio <- t(IENCEratio)
IENCEratio <- data.frame(cluster=row.names(IENCEratio),IENCEratio)
IENCEratio$cluster <- gsub("X","C",IENCEratio$cluster)
IENCEratio$cluster <- factor(IENCEratio$cluster,levels=as.character(IENCEratio$cluster))

sampleCompWideNCIEInd <- sampleCompWide/apply(sampleCompWide,1,sum)
sampleCompWideNCIEInd <- data.frame(type = c("IE","IE","NCE","IE","NCE","NCE","NCE","IE","IE"), sampleCompWideNCIEInd)
sampleCompWideNCIEInd <- gather(sampleCompWideNCIEInd, cluster, proportion, X0:X12, factor_key=TRUE)
sampleCompWideNCIEInd$cluster <- gsub("X","C",sampleCompWideNCIEInd$cluster)
sampleCompWideNCIEInd$cluster <- factor(sampleCompWideNCIEInd$cluster,levels=as.character(IENCEratio$cluster))
  
plot <- ggplot(IENCEratio) +
  geom_bar(aes(x = cluster,y = log(IE,2),fill=cluster), stat='identity') +
  geom_hline(yintercept = 0,size=0.25) +
  labs(title = "Combined Composition", x="Cluster", y="Log2(IE:NCE)") +
  theme_base() +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=6),text=element_text(size=6),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_rect(colour='black'))
ggsave(paste("nc_combined_2/sampleComp_combined_IEvsNCE_by_subgroup",".pdf",sep=''), plot, width = 2.5, height = 2, useDingbats = FALSE, limitsize=FALSE)


plot <- ggplot(sampleCompWideNCIEInd) +
  geom_boxplot(aes(x = cluster, y = proportion, fill = type),outlier.size=0.25,size=0.25) +
  labs(title = "Sample Composition", x="Cluster", y="Proportion of sample") +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_rect(colour='black'))
ggsave(paste("nc_combined_2/sampleComp_proportion_by_subgroup",".pdf",sep=''), plot, width = 3, height = 2, useDingbats = FALSE, limitsize=FALSE)

c = "C12"
t.test(sampleCompWideNCIEInd[sampleCompWideNCIEInd$type == "IE" & sampleCompWideNCIEInd$cluster == c,3],
       sampleCompWideNCIEInd[sampleCompWideNCIEInd$type == "NCE" & sampleCompWideNCIEInd$cluster == c,3])

# diversity index
sampleCompDiversity <- data.frame(type = c("IE","IE","NCE","IE","NCE","NCE","NCE","IE","IE"), 
                                  simpson = apply(sampleCompWide,1,function(x) x = diversity(x, "simpson")))

plot <- ggplot(sampleCompDiversity) +
  geom_boxplot(aes(x="",y = simpson, fill = type),outlier.size=0.25,size=0.25) +
  labs(title = "", x="", y="Simpson Diversity Index") +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_rect(colour='black'))
ggsave(paste("nc_combined_2/simpson_diversity_subgroup",".pdf",sep=''), plot, width = 2, height = 2, useDingbats = FALSE, limitsize=FALSE)


###################################################
# composite signatures for rat injury model genes #
k27Cut <- as.character(read.table('rat_H3K27ac_InjuryDB_closest_genes_specific_human.txt',header=FALSE)$V1)
k27Cut <- k27Cut[!is.na(k27Cut)]

k27CutExpr <- intersect(row.names(GetAssayData(dat.all.filtered,slot = 'scale.data')),k27Cut)
k27CutExprTop50 <- intersect(k27CutExpr,as.character(top50$gene))

# gather nerve injury score for clusters
k27CutExprDat <- data.frame(embeddings, k27Cut = apply(GetAssayData(dat.all.filtered,slot = 'scale.data')[k27CutExpr,row.names(embeddings)],2,mean),
                                        k27CutTop50 = apply(GetAssayData(dat.all.filtered,slot = 'scale.data')[k27CutExprTop50,row.names(embeddings)],2,mean))

k27CutExprDat <- k27CutExprDat[order(k27CutExprDat$k27Cut),]

# plot the mean nerve injury signatures
plot <- ggplot(k27CutExprDat) +
  geom_point(aes(x = UMAP_1, y = UMAP_2, color = k27Cut), size=0.25) +
  scale_color_gradientn(colors = rev(brewer.pal(8, 'RdBu'))) +
  labs(title = "Rat Injury K27 signature", x="UMAP 1", y="UMAP 2") +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=0,hjust=0,vjust=1),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_blank(),panel.background=element_blank(),axis.line = element_line(size=0.25))
ggsave("nc_combined_2/vs_nc_combined_2_rat_injury_k27.png", plot, width = 5, height = 4)


k27CutExprDat <- k27CutExprDat[order(k27CutExprDat$k27Cut),]

# plot the mean nerve injury signatures of top 50 marker genes
k27CutExprDat <- k27CutExprDat[order(k27CutExprDat$k27CutTop50),]

plot <- ggplot(k27CutExprDat) +
  geom_point(aes(x = UMAP_1, y = UMAP_2, color = k27CutTop50), size=0.25) +
  scale_color_gradientn(colors = rev(brewer.pal(8, 'RdBu'))) +
  labs(title = "Rat Injury K27 signature - overlap with top 50 markers", x="UMAP 1", y="UMAP 2") +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=90,hjust=0,vjust=1),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_blank(),panel.background=element_blank(),axis.line = element_line(size=0.25))
ggsave("nc_combined_2/vs_nc_combined_2_rat_injury_k27_top50.png", plot, width = 5, height = 4)

write.table(k27CutExprDat,'nc_combined_2/vs_nc_combined_2_UMAP_res0.3_rat_injury_embeddings.txt',sep='\t',quote=FALSE,row.names=TRUE,col.names=NA)

# SCPubr featureplot
modulescored <- AddModuleScore(dat.all.filtered, list(k27CutExpr),name = "Nerve_Injury")
SCpubr::do_FeaturePlot(sample = modulescored, features = "Nerve_Injury1",order=TRUE,enforce_symmetry = TRUE)
ggsave(paste("nc_combined_2/featureplots/vs_nc_combined_2_umap_feature_","Nerve_Injury","_new.pdf",sep=''), width = 5, height = 6)


###################################
## RT compared to non RT samples ##

# data filter
dat.RTcomp <- dat.all.filtered

# subset to those samples of interest and normalize the data
Idents(dat.RTcomp) <- dat.RTcomp$orig.ident
dat.RTcomp <- subset(x = dat.RTcomp, idents = c("S11", "S12","S29","VSF02","VSF03")) # all from separate batches 

dat.RTcomp$RT <- dat.RTcomp$orig.ident
dat.RTcomp$RT <- gsub("^S11$","RT",dat.RTcomp$RT)
dat.RTcomp$RT <- gsub("^S12$","RT",dat.RTcomp$RT)
dat.RTcomp$RT <- gsub("^S29$","noRT",dat.RTcomp$RT)
dat.RTcomp$RT <- gsub("^VSF02$","noRT",dat.RTcomp$RT)
dat.RTcomp$RT <- gsub("^VSF03$","noRT",dat.RTcomp$RT)

# SC transform - all 3 samples for comparison are from 3 batches
dat.RTcomp <- SCTransform(object = dat.RTcomp, vars.to.regress = c("nCount_RNA", "percent.mito","orig.ident"), verbose = FALSE)

# linear dimensionality reduction
dat.RTcomp <- RunPCA(object = dat.RTcomp, verbose = FALSE)
DimPlot(object = dat.RTcomp, reduction = "pca",pt.size = 0.25,group.by = 'orig.ident')
ElbowPlot(object = dat.RTcomp,ndims = 30) 

# run harmony through seurat
dat.RTcomp <- RunHarmony(dat.RTcomp, c("source","batch"),max.iter.harmony = 20, plot_convergence = TRUE,theta = 1)

# Run UMAP
mindist=0.3
dat.RTcomp <- RunUMAP(dat.RTcomp, reduction = "harmony",dims = 1:12, min.dist = mindist)
plot <- DimPlot(object = dat.RTcomp, reduction = "umap", group.by = "orig.ident", pt.size=0.3, do.label=T)
ggsave(paste("nc_combined_2/rt_comp/vs_snRNA_RTnoRT_harmony_umap_mindist",mindist,".png",collapse="",sep=""), plot, width =6, height = 5)

# find clusters using harmony embeddings
dat.RTcomp <- FindNeighbors(object = dat.RTcomp, dims=1:12,reduction = "harmony")
res=0.2
dat.RTcomp <- FindClusters(object = dat.RTcomp, reduction = "harmony", resolution = res)

# plot clusters in UMAP space 
plot <- DimPlot(object = dat.RTcomp, reduction = "umap", pt.size=0.3, do.label=T, group.by = paste("SCT_snn_res.",res,collapse="",sep=""))
ggsave(paste("nc_combined_2/rt_comp/vs_snRNA_RTnoRT_harmony_umap_clusters_mindist",mindist,"_res",res,".png",collapse="",sep=""), plot, width =5, height = 4)

# identify cluster markers
vsRTnoRT.markers <- FindAllMarkers(object = dat.RTcomp, only.pos = TRUE, min.pct = 0.1, thresh.use = 0.1)
vsRTnoRT.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)

# marker heatmap
top10 <- vsRTnoRT.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
write.table(top10,paste("nc_combined_2/rt_comp/top10markers_snRNA_RTnoRT_harmony_res",res,".txt",sep=''),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)

# top 50
top50 <- vsRTnoRT.markers %>% group_by(cluster) %>% top_n(50, avg_logFC)
write.table(top50,paste("nc_combined_2/rt_comp/top50markers_snRNA_RTnoRT_harmony_res",res,".txt",sep=''),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)

plot <- DoHeatmap(object = dat.RTcomp, features = intersect(top10$gene,row.names(dat.RTcomp@assays$SCT@scale.data)), slot = 'scale.data',draw.lines = FALSE)
ggsave(paste("nc_combined_2/rt_comp/vs_snRNA_RTnoRT_harmony_heatmap_top10marker_scaled_res",res,".pdf",sep=''), plot, width = 9, height = 12, useDingbats = FALSE, limitsize=FALSE)
ggsave(paste("nc_combined_2/rt_comp/vs_snRNA_RTnoRT_harmony_heatmap_top10marker_scaled_res",res,".png",sep=''), plot, width = 9, height = 12,limitsize = FALSE)

# ordered featureplots
goi <- c("APOE","APOC1","SORCS2")
for (i in 1:length(goi)){
  fp <- FeaturePlot(object = dat.RTcomp, reduction = "umap", features = goi[i], pt.size = 0.1)
  fp$data <- fp$data[order(fp$data[,3]),]
  ggsave(paste("nc_combined_2/rt_comp/vs_nc_combined_RTnoRT_umap_feature_",goi[i],".pdf",sep=''), fp, width = 3, height = 2.5, useDingbats = FALSE, limitsize=FALSE)
}

##################################################
## AddModuleScore requests from Nature reviewer ##
# obtain gene sets that make sense
GO_CytokineSig <- as.character(read.table('~/common/genesets/GO_CYTOKINE_MEDIATED_SIGNALING_PATHWAY.txt',sep='\t',skip=2,header=FALSE)$V1)
GO_EnsheathmentNeuron <- as.character(read.table('~/common/genesets/GO_ENSHEATHMENT_OF_NEURONS.txt',sep='\t',skip=2,header=FALSE)$V1)
GO_BloodVessel <- as.character(read.table('~/common/genesets/GO_BLOOD_VESSEL_MORPHOGENESIS.txt',sep='\t',skip=2,header=FALSE)$V1)
LEE_DIFFERENTIATING_T_LYMPHOCYTE <- as.character(read.table('~/common/genesets/LEE_DIFFERENTIATING_T_LYMPHOCYTE.txt',sep='\t',skip=2,header=FALSE)$V1)
GO_NEUROGENESIS <- as.character(read.table('~/common/genesets/GO_NEUROGENESIS.txt',sep='\t',skip=2,header=FALSE)$V1)
GO_NEURON_PROJECTION <- as.character(read.table('~/common/genesets/GO_NEURON_PROJECTION.txt',sep='\t',skip=2,header=FALSE)$V1)  
Zhao_Endothelial <- as.character(read.table('~/common/genesets/zhao_2018_EC_markers.txt',sep='\t',skip=2,header=FALSE)$V1)  # Zhao et al. 2018; Cancer Res
GO_INFLAMMATORY_RESPONSE <- as.character(read.table('~/common/genesets/GO_INFLAMMATORY_RESPONSE.txt',sep='\t',skip=2,header=FALSE)$V1)  
GO_COMPLEMENT_ACTIVATION <- as.character(read.table('~/common/genesets/GO_COMPLEMENT_ACTIVATION.txt',sep='\t',skip=2,header=FALSE)$V1)  
BOQUEST_STEM_CELL_UP <- as.character(read.table('~/common/genesets/BOQUEST_STEM_CELL_UP.txt',sep='\t',skip=2,header=FALSE)$V1)  
GO_GLIOGENESIS <- as.character(read.table('~/common/genesets/GO_GLIOGENESIS.txt',sep='\t',skip=2,header=FALSE)$V1)  
EBAUER_MYOGENIC_TARGETS_OF_PAX3_FOXO1_FUSION <- as.character(read.table('~/common/genesets/EBAUER_MYOGENIC_TARGETS_OF_PAX3_FOXO1_FUSION.txt',sep='\t',skip=2,header=FALSE)$V1)  
microglia_manual <- c("P2Y12","TMEM119","C1QB","C1QC") # manual curation

# addmodulescore function
aset = microglia_manual
asetName = "microglia_markers"
bset = 12

# modulescored <- AddModuleScore(dat.all.filtered, list(intersect(aset,as.character(data.frame(top50)[top50$cluster == bset,"gene"]))),name = paste(bset,asetName,sep = "_"))
# fp <- FeaturePlot(object = modulescored, reduction = "umap", features = paste(paste(bset,asetName,sep = "_"),1,sep=""), pt.size = 0.1)
# fp$data <- fp$data[order(fp$data[,3]),]
# ggsave(paste("nc_combined_2/modulescorePlots/vs_nc_combined_2_umap_ModuleScore_",paste(bset,asetName,sep = "_"),".pdf",sep=''), fp, width = 3, height = 2.5, useDingbats = FALSE, limitsize=FALSE)

# New featureplots 11/7/22
modulescored <- AddModuleScore(dat.all.filtered, list(intersect(aset,as.character(data.frame(top50)[top50$cluster == bset,"gene"]))),name = paste(bset,asetName,sep = "_"))
SCpubr::do_FeaturePlot(sample = modulescored, features = paste("X",paste(bset,asetName,sep = "_"),1,sep=""),order=TRUE)
ggsave(paste("nc_combined_2/modulescorePlots/vs_nc_combined_2_umap_ModuleScore_",paste(bset,asetName,sep = "_"),"new.pdf",sep=''), width = 5, height = 6)

#Manual score
modulescored <- AddModuleScore(dat.all.filtered, list(aset),name = asetName)
fp <- FeaturePlot(object = modulescored, reduction = "umap", features = paste(asetName,1,sep=""), pt.size = 0.1)
fp$data <- fp$data[order(fp$data[,3]),]
ggsave(paste("nc_combined_2/modulescorePlots/vs_nc_combined_2_umap_ModuleScore_",asetName,".pdf",sep=''), fp, width = 3, height = 2.5, useDingbats = FALSE, limitsize=FALSE)

# New featureplots for manual 11/7/22
modulescored <- AddModuleScore(dat.all.filtered, list(aset),name = asetName)
SCpubr::do_FeaturePlot(sample = modulescored, features = paste(asetName,1,sep=""),order=TRUE)
ggsave(paste("nc_combined_2/modulescorePlots/vs_nc_combined_2_umap_ModuleScore_",asetName,"new.pdf",sep=''), width = 5, height = 6)

# Cluster inter-replicate variability.
sampleCompWideAnnot <- sampleCompWide
sampleCompWideAnnot <- sampleCompWideAnnot/apply(sampleCompWideAnnot,1,sum)
sampleCompWideAnnot <- round(sampleCompWideAnnot,5)
sampleCompWideAnnot <- data.frame(sampleCompWideAnnot,subgroup = c("I","I","N","I","N","N","N","I","I"))

pVals <- data.frame(cluster = colnames(sampleCompWideAnnot)[1:13], pval= rep(0,13))

for(i in 1:13){
  clust = colnames(sampleCompWideAnnot)[i]
  dist1 = sampleCompWideAnnot[sampleCompWideAnnot$subgroup == "I",clust]
  dist2 = sampleCompWideAnnot[sampleCompWideAnnot$subgroup == "N",clust]
  pVals[i,2] <- t.test(dist1,dist2,var.equal = TRUE, paired = FALSE)[['p.value']]
}

# Violin plot for Gli3 by request 4/14/20
plot <- VlnPlot(object = dat.all.filtered, features = "GLI3", pt.size = 0.1)
ggsave(paste("nc_combined_2/violinplots/vs_nc_combined_2_violin_","GLI3",".pdf",sep=''), plot, width = 3, height = 2.5, useDingbats = FALSE, limitsize=FALSE)

## Microglia subtype analysis from Chang Kim in Tom Nowakowski's lab 8/13/20
kimRDS <- readRDS("swan.rds")
kimEmbeddings <- data.frame(embeddings)
row.names(kimEmbeddings) <- paste(kimEmbeddings$sample,row.names(kimEmbeddings),sep = "_")  
row.names(kimEmbeddings) <- paste(row.names(kimEmbeddings),"-1",sep="")
kimEmbeddings <- data.frame(kimEmbeddings, kimMicroglia=as.character((kimRDS$clusters[row.names(kimEmbeddings)])))

plot <- ggplot(kimEmbeddings) +
  geom_point(aes(x = UMAP_1, y = UMAP_2, color = kimMicroglia), size=0.25) +
  labs(title = "Microglia subtypes", x="UMAP 1", y="UMAP 2") +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=0,hjust=0,vjust=1),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_blank(),panel.background=element_blank(),axis.line = element_line(size=0.25))
ggsave("nc_combined_2/vs_nc_combined_2_Kim_microglia.png", plot, width = 5, height = 4)

# Continue to look at candidate genes
SCpubr::do_FeaturePlot(sample = dat.all.filtered, features = "PTPRG", max.cutoff = 1.5,order=TRUE)
ggsave("nc_combined_2/vs_nc_combined_2_PTPRG.pdf", width = 5, height = 6)

# Violin plots on request 
VlnPlot(object = dat.all.filtered, features = "PTPRZ1", pt.size = 0.1)
SCpubr::do_ViolinPlot(sample = dat.all.filtered, features = "PTPRZ1")
SCpubr::do_RidgePlot(sample = dat.all.filtered,feature = "PTPRZ1")

# Load TCF3 motifs from snARC-seq results
SETDB1_perturb_at_TCF3 <- as.character(read.table('/raleighlab/data1/liuj/schwannoma/snarc-seq/analysis/snarc_rna_2/SETDB1_perturb_at_TCF3.txt',sep='\t',header=FALSE)$V1)  

modulescored <- AddModuleScore(dat.all.filtered, list(intersect(SETDB1_perturb_at_TCF3,row.names(dat.all.filtered))),name = "SETDB1_perturb_at_TCF3")
SCpubr::do_FeaturePlot(sample = modulescored, features = "SETDB1_perturb_at_TCF31",order=TRUE)
ggsave("nc_combined_2/vs_nc_combined_2_featureplot_SETDB1_perturb_at_TCF3.pdf", width = 5, height = 6)

SCpubr::do_ViolinPlot(sample = modulescored, features = "SETDB1_perturb_at_TCF31")
ggsave("nc_combined_2/vs_nc_combined_2_violin_SETDB1_perturb_at_TCF3.pdf", width = 6, height = 3)
