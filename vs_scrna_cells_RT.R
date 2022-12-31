## Analysis of HEI193 scRNA-seq of irradiated cells in triplicates
library(Seurat)
library(ggplot2)
library(ggthemes)

setwd('/raleighlab/data1/liuj/schwannoma/scrna_cells/analysis')

# read in data
smpid <- dir(file.path("/media/data2/liuj/vestibular/scrna_cells/cellranger_out/"))
dataGene.list <- list()

for (i in 1:length(smpid)){
  df <- Read10X(data.dir = paste("/media/data2/liuj/vestibular/scrna_cells/cellranger_out/",smpid[i],"/outs/filtered_feature_bc_matrix/",sep=''))
  dataGene.list[[i]] <- CreateSeuratObject(df, project = smpid[i], min.features = 100,min.cells = 5)
}

names(dataGene.list) <- smpid

# merge samples 
dat.all <- merge(x = dataGene.list[[smpid[1]]], y = dataGene.list[smpid[2:length(smpid)]])

# mitochondrial analysis
mito.features <- grep(pattern = "^MT-", x = rownames(x = dat.all), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = dat.all, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = dat.all, slot = 'counts'))
dat.all[['percent.mito']] <- percent.mito

plot <- VlnPlot(object = dat.all, features= c("nFeature_RNA", "nCount_RNA", "percent.mito"),pt.size = 0,group.by = 'orig.ident')
ggsave("plots/vs_scrna_cells_qc.png", plot, width = 10, height = 6)

# remove low quality cells 
dat.all <- subset(dat.all, subset = nCount_RNA > 2000 & percent.mito < 0.20)
dat.all <- SCTransform(dat.all, vars.to.regress = c("percent.mito"), verbose = FALSE)
dat.all <- RunPCA(object = dat.all, verbose = FALSE)

ElbowPlot(object = dat.all)

# Run UMAP embeddings
mindist=0.2
dat.all <- RunUMAP(dat.all, reduction = "pca",dims = 1:30, min.dist = mindist)
plot <- DimPlot(object = dat.all, reduction = "umap", group.by = "orig.ident", pt.size=0.3)
ggsave("plots/umap_sct.png", plot, width =8, height = 6)

plot <- DimPlot(object = dat.all, reduction = "umap", split.by  = "orig.ident", pt.size=0.3)
ggsave("plots/umap_sct_split.png", plot, width =24, height = 6)

# plot by condition
dat.all$cond <- dat.all$orig.ident
dat.all$cond <- gsub(pattern = "_1$","",dat.all$cond)
dat.all$cond <- gsub(pattern = "_2$","",dat.all$cond)
dat.all$cond <- gsub(pattern = "_3$","",dat.all$cond)

plot <- DimPlot(object = dat.all, reduction = "umap", group.by  = "cond", pt.size=0.3)
ggsave("plots/umap_sct_cond.png", plot, width =8, height = 6)

# find clusters, what is the heterogeneity due to?
dat.all <- FindNeighbors(object = dat.all, dims = 1:30) 

res=0.5
dat.all <- FindClusters(object = dat.all, reduction = "umap", resolution = res)

plot <- DimPlot(object = dat.all, reduction = "umap", pt.size=0.3, group.by = "SCT_snn_res.0.3")
ggsave(paste("plots/umap_sct_clusters_res",res,".png",sep=""), plot, width =6.5, height = 6)

plot <- DimPlot(object = dat.all, reduction = "umap", pt.size=0.3, group.by = "SCT_snn_res.0.3", label = TRUE)
ggsave(paste("plots/umap_sct_clusters_res",res,".png",sep=""), plot, width =6.5, height = 6)

plot <- DimPlot(object = dat.all, reduction = "umap", split.by  = "cond", pt.size=0.3, group.by = "SCT_snn_res.0.3")
ggsave(paste("plots/umap_sct_cond_split_res",res,".png",sep=""), plot, width =16, height = 6)

plot <- DimPlot(object = dat.all, reduction = "umap", split.by  = "orig.ident", pt.size=0.3, group.by = "SCT_snn_res.0.3")
ggsave(paste("plots/umap_sct_orig_ident_split_res",res,".png",sep=""), plot, width =24, height = 6)

# identify cluster markers
markers <- FindAllMarkers(object = dat.all, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25, assay="SCT",slot='data',test.use='wilcox',pseudocount.use = 1)
markers <- FindAllMarkers(object = dat.all, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

markers %>% group_by(cluster) %>% top_n(2, avg_log2FC)

# marker heatmap
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)

top50 <- markers %>% group_by(cluster) %>% top_n(50, avg_log2FC)
write.table(top50,paste("plots/top50markers_res",res,".txt",sep=''),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)

plot <- DoHeatmap(subset(dat.all, downsample = 200), features = top10$gene, size = 3) + NoLegend()
ggsave(filename = paste("plots/vs_HEI193_RT_scrna_top10_heatmap_15clusters_res",res,".pdf", sep=''),plot, width = 8, height = 18)

#df <- plot$data[,c("Feature","Expression","Identity")]
#df <- reshape(df, idvar = "Feature", timevar = "Identity", direction = "wide")
write.table(plot$data,paste("plots/vs_HEI193_RT_scrna_top10_heatmap_15clusters_res",res,".txt",sep=''),sep='\t',quote=FALSE,row.names=FALSE)

# add cell cycle score
dat.all <- CellCycleScoring(dat.all, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)

plot <- VlnPlot(dat.all, group.by = 'SCT_snn_res.0.5', features=c('S.Score','G2M.Score'),pt.size = 0)
ggsave("plots/vs_scrna_violin_cellcycle.png", plot, width = 4, height = 3)

plot <- FeaturePlot(dat.all, features = c('S.Score','G2M.Score'), split.by = "cond",cols = c("grey", "purple"),pt.size = 0.25,order=TRUE,min.cutoff = 'q5', max.cutoff = 'q80')
ggsave("featureplots/vs_scrna_featureplot_cellcycle.png", plot, width = 9, height = 7)

# featureplots to start 
goi = grep(pattern = "^APO", x = rownames(x = dat.all), value = TRUE)
goi = grep("APOB",goi,invert = TRUE, value=TRUE)

plot <- FeaturePlot(dat.all, features = goi, split.by = "cond",cols = c("grey", "purple"),pt.size = 0.25,order=TRUE)
ggsave("featureplots/vs_scrna_featureplot_APO.png", plot, width = 12, height = 38)

plot <- VlnPlot(dat.all, group.by = 'SCT_snn_res.0.5', features=goi,pt.size = 0.5)
ggsave("plots/vs_scrna_violineplot_APO_cluster.png", plot, width = 12, height = 20)

# look at screen hits 
goi <- c('KDM3B', 'ELP5', 'YEATS4', 'KANSL3', 'GPS2', 'TADA1', 'ENY2', 'MRGBP', 'KAT7', 'ELP2', 'SUDS3', 'KANSL2', 'KDM1A', 'MSL3', 'KAT8', 'SETDB1', 'DMAP1', 'BRPF1', 'RPS2', 'SUPT20H', 'KDM5C', 'TRRAP', 'ELP4', 'EHMT1', 'ZZZ3', 'KDM4A', 'PRMT5', 'ING5', 'TAF6L', 'TADA2B', 'YEATS2', 'CHD4', 'SUV39H1', 'EZH2', 'ASH2L', 'ELP6', 'CREBBP', 'SMARCA2', 'SUPT7L', 'HDAC8', 'EP400', 'CDK4', 'TAF10', 'HAT1', 'HMG20B', 'OGT', 'EPC1', 'ING3', 'SUZ12', 'EED', 'NF2', 'SMARCB1', 'LZTR1')

plot <- FeaturePlot(dat.all, features = goi, split.by = "cond",cols = c("grey", "purple"),pt.size = 0.25,order=TRUE)
ggsave("featureplots/vs_scrna_featureplot_screen_hits.png", plot, width = 8, height = 100, limitsize = FALSE)

# get expressed genes list for motif subsetting
df <- GetAssayData(dat.all, slot="counts")
expressedGenes <- row.names(df[apply(df,1,function(x) max(x) > 10),])
write.table(expressedGenes,'expressedGenes_in_HEI193_max10.txt',sep='\t',quote=FALSE,row.names=FALSE)

# export cluster composition
write.table(clusterComp02,'cluster_composition_HEI193_RT_res05.txt',sep='\t',quote=FALSE,row.names=FALSE)

# export embeddings for further exploration
embeddings <- Embeddings(dat.all, reduction="umap")
embeddings <- data.frame(embeddings, cluster=as.character((dat.all$SCT_snn_res.0.5[row.names(embeddings)])),
                         condition=as.character((dat.all$cond[row.names(embeddings)])),
                         sample=as.character((dat.all$orig.ident[row.names(embeddings)])))
write.table(embeddings,'HEI193_RT_UMAP_res0.5_embeddings.txt',sep='\t',quote=FALSE,row.names=TRUE,col.names=NA)

