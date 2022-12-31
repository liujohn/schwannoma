## Analysis of in vitro HEI-193 perturb-seq 
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(ggthemes)
library(plyr)
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(DESeq2)
library(pheatmap)
library(tuple)
library(pals)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)
setwd('/raleighlab/data1/liuj/schwannoma/perturb-seq/analysis/perturb_2/')


# load all vs in vitro data
smpid <- dir(file.path("/raleighlab/data1/liuj/schwannoma/perturb-seq/cellranger_out/perturb_1/"))

dataGene.list <- list()
sgR.list <- list()

for (i in 1:length(smpid)){
  df <- Read10X(data.dir = paste("/raleighlab/data1/liuj/schwannoma/perturb-seq/cellranger_out/perturb_1/",smpid[i],"/outs/filtered_feature_bc_matrix/",sep=''))
  dataGene.list[[i]] <- CreateSeuratObject(counts=df$`Gene Expression`, project = smpid[i])
  sgR.list[[i]] <- CreateSeuratObject(counts = df$`CRISPR Guide Capture`, project = smpid[i])
  dataGene.list[[i]]$orig.ident <- smpid[i]
}

names(dataGene.list) <- smpid
names(sgR.list) <- smpid

dat.all <- merge(x = dataGene.list[[1]], y = dataGene.list[2:length(dataGene.list)])

# quality metrics for samples divided by ID 
mito.features <- grep(pattern = "^MT-", x = rownames(x = dat.all), value = TRUE, ignore.case = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = dat.all, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = dat.all, slot = 'counts'))
dat.all[['percent.mt']] <- percent.mito

VlnPlot(dat.all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
ggsave(filename = "qc_violin_sub.pdf",width = 16, height=10)

# subset data and normalize
dat.sub <- subset(dat.all, subset = nFeature_RNA > 200)
dat.sub <- SCTransform(dat.sub, verbose = TRUE)
dat.sub <- RunPCA(dat.sub, verbose = TRUE)

# PCA dimensions
ElbowPlot(object = dat.sub,ndims = 50) 

# assign annotations to samples
dat.all$RT <- dat.all$orig.ident
dat.all$RT <- gsub("VS_perturb_1","0Gy",dat.all$RT)
dat.all$RT <- gsub("VS_perturb_2","0Gy",dat.all$RT)
dat.all$RT <- gsub("VS_perturb_3","9Gy",dat.all$RT)
dat.all$RT <- gsub("VS_perturb_4","12_5Gy",dat.all$RT)

dat.sub$RT <- dat.sub$orig.ident
dat.sub$RT <- gsub("VS_perturb_1","0Gy",dat.sub$RT)
dat.sub$RT <- gsub("VS_perturb_2","0Gy",dat.sub$RT)
dat.sub$RT <- gsub("VS_perturb_3","9Gy",dat.sub$RT)
dat.sub$RT <- gsub("VS_perturb_4","12_5Gy",dat.sub$RT)

# Run UMAP embeddings
mindist=0.7
dat.sub <- RunUMAP(dat.sub, reduction = "pca",dims = 1:30, min.dist = mindist)
DimPlot(object = dat.sub, group.by = 'orig.ident', reduction = "umap", pt.size=0.3,raster=FALSE)
ggsave("pheno_umap_sct_sub.png", width =14, height = 9)

DimPlot(object = dat.sub, group.by = 'RT', reduction = "umap", pt.size=0.3,raster=FALSE)
ggsave("pheno_umap_sct_sub_RT.png", width =12, height = 9)

# find clusters, what is the heterogeneity due to?
dat.sub <- FindNeighbors(object = dat.sub, dims = 1:30) 

res=0.4
dat.sub <- FindClusters(object = dat.sub, reduction = "umap", resolution = res)

# plot clusters in UMAP space 
DimPlot(object = dat.sub, reduction = "umap", pt.size=0.3, group.by = paste("SCT_snn_res.",res,collapse="",sep=""),label = TRUE,raster = FALSE)
ggsave(paste("pheno_sct_sub_umap_mindist_",mindist,"_res",res,".png",collapse="",sep=""), width =12, height = 9)


##########################
# Generate sgRNA calls 
# Fill in sgRNA metadata 
sgR.all <- merge(x = sgR.list[[1]], y = sgR.list[2:length(sgR.list)])

sgR.sub <- subset(sgR.all, nFeature_RNA > 0)
sgR.sub <- GetAssayData(sgR.sub, slot="counts")
sgR.sub <- data.frame(sgR.sub)
sgR.sub <- sgR.sub[apply(sgR.sub,1,sum) > 0,]
colnames(sgR.sub) <- gsub("\\.","-",colnames(sgR.sub))

seqMetrics.list <- list()
sgRNATbl.list <- list()

for (i in 1:length(smpid)){
  seqMetrics.list[[i]] <- read.table(paste("/raleighlab/data1/liuj/schwannoma/perturb-seq/cellranger_out/perturb_1/",smpid[i],"/outs/metrics_summary.csv",sep=''),sep=',',header=TRUE)
  sgRNA_tags <- read.table(paste("/raleighlab/data1/liuj/schwannoma/perturb-seq/cellranger_out/perturb_1/",smpid[i],"/outs/crispr_analysis/protospacer_calls_per_cell.csv",sep=''),sep=',',header=TRUE)
  sgRNA_tags$cell_barcode <- gsub("-1","",sgRNA_tags$cell_barcode)
  row.names(sgRNA_tags) <- as.character(sgRNA_tags$cell_barcode)
  sgRNATbl.list[[i]] <- sgRNA_tags
}


# histogram of sgRNA dectections
names(sgRNATbl.list) <- smpid
sgRNATbl <- ldply (sgRNATbl.list, data.frame)

plot <- ggplot(sgRNATbl) + 
  geom_histogram(aes(x=num_features, fill=.id),color=NA,show.legend = FALSE,size=0.25,binwidth=1) +
  scale_x_continuous(limits = c(0,10),breaks = seq(0,10,1),expand=c(0.01,0.01)) +
  scale_y_continuous(expand=c(0.01,0.01), limits = c(0,3500)) +
  labs(title = "Number of sgRNAs detected", x="sgRNAs", y="cells") +
  facet_wrap(~.id, scales="free", ncol = 4) +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=0,hjust=0.5,vjust=1),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_rect(colour='black',fill=NA, size=0.25),panel.background=element_blank(),axis.line = element_line(size=0.25))
ggsave("pdx_perturb_1_unique_sgRNAs_hist_sub.png", plot, width =6, height = 2)


# number of concordant sgRNAs 
sgRNATblConcordant <- sgRNATbl[sgRNATbl$num_features == 1,]
sgRNATblConcordant$Gene <- sapply(strsplit(sgRNATblConcordant$feature_call,"_"), `[`, 1)

# heatmap visualization of concordant sgRNA detections 
df <- melt(table(sgRNATblConcordant[,c(".id","feature_call")]))

plot <- ggplot(df, aes(x = feature_call,.id)) +
  geom_tile(aes(fill = value)) +
  geom_text(aes(label = value),size=1) +
  scale_fill_gradient(low = "white", high = "steelblue",limit=c(0,40),na.value="steelblue") +
  labs(title = "Number of cells with single sgRNAs detected", x="", y="") +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_blank(),panel.background=element_blank(),axis.line = element_line(size=0.25))
ggsave("pdx_perturb_vs_sgRNAs_heatmap_sub.png", plot, width =9, height = 1.8)


plot <- ggplot(df) + 
  geom_histogram(aes(x=value,fill=.id),color=NA,show.legend = FALSE,size=0.25,binwidth=10) +
  scale_x_continuous(limits = c(-20,100), expand=c(0.01,0.01)) +
  #scale_y_continuous(limits = c(0,51), expand=c(0.01,0.01)) +
  labs(title = "N cells per sgRNA pair detected", x="Cells", y="Perturbations (single sgRNA)") +
  facet_wrap(~.id,scales="fixed", ncol = 4) +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=0,hjust=0.5,vjust=1),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_rect(colour='black',fill=NA, size=0.25),panel.background=element_blank(),axis.line = element_line(size=0.25))
ggsave("pdx_perturb_cells_per_sgRNA_hist_sub.png", plot, width =6, height = 2)

# get UMAP coordinates for sgRNA feature plots - concordant cells only
subEmbeddings <- data.frame(Embeddings(object = dat.sub, reduction = "umap"))
subEmbeddings$orig.ident <- dat.sub$orig.ident[row.names(subEmbeddings)]
subEmbeddings$RT <- dat.sub$RT[row.names(subEmbeddings)]
write.table(subEmbeddings,'vs_perturb_reanalysis_102022_embeddings.txt',sep='\t',quote=FALSE,row.names=TRUE,col.names=NA)

# remake sgRNA feature plots using concordant cells. 
for(i in 1:length(smpid)){
  df = subEmbeddings[subEmbeddings$orig.ident == smpid[[i]],]
  row.names(df) = sapply(strsplit(row.names(df),"-"), `[`, 1)
  whitelistCells = sgRNATblConcordant[sgRNATblConcordant$.id  == smpid[[i]],"cell_barcode"]
  dfWhiteListed = df[whitelistCells,]
  
  Npos = nrow(dfWhiteListed)
  Ntot = nrow(df)
  
  plot <- ggplot() +
    geom_point(data= df, aes(x = UMAP_1, y = UMAP_2), color = "gray50", size=0.25) +
    geom_point(data= dfWhiteListed, aes(x = UMAP_1, y = UMAP_2), color='red', size=0.25) +
    scale_color_gradientn(colors = brewer.pal(8, 'Reds'),) +
    #scale_x_continuous(limits = c(-16,16)) +
    #scale_y_continuous(limits = c(-16,16)) +
    labs(title = paste(smpid[[i]],"\n",Npos,"sgRNA+ cells /",Ntot,"total cells",round(Npos/Ntot*100,3),"% Pass QC"), x="UMAP 1", y="UMAP 2") +
    theme_base() +
    theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=0,hjust=0,vjust=1),
          axis.ticks=element_line(colour='black',size=0.25),panel.border=element_blank(),panel.background=element_blank(),axis.line = element_line(size=0.25))
  ggsave(paste("vs_sgRNA_reads_binary_concordant_",smpid[[i]],".png",sep=""), plot, width =5, height = 5)
}

# combined sgRNA plot 
df = subEmbeddings
row.names(df) = sapply(strsplit(row.names(df),"-"), `[`, 1)
whitelistCells = sgRNATblConcordant[,"cell_barcode"]
dfWhiteListed = df[whitelistCells,]

Npos = nrow(dfWhiteListed)
Ntot = nrow(df)

plot <- ggplot() +
  geom_point(data= df, aes(x = UMAP_1, y = UMAP_2), color = "gray50", size=0.25) +
  geom_point(data= dfWhiteListed, aes(x = UMAP_1, y = UMAP_2), color='red', size=0.25) +
  scale_color_gradientn(colors = brewer.pal(8, 'Reds'),) +
  scale_x_continuous(limits = c(-16,16)) +
  scale_y_continuous(limits = c(-16.8,16)) +
  labs(title = paste("vs","\n",Npos,"sgRNA+ cells /",Ntot,"total cells",round(Npos/Ntot*100,3),"% Pass QC"), x="UMAP 1", y="UMAP 2") +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=0,hjust=0,vjust=1),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_blank(),panel.background=element_blank(),axis.line = element_line(size=0.25))
ggsave(paste("vs_sgRNA_reads_binary_concordant_","combined",".png",sep=""), plot, width =6, height = 5)


#########################
# Get knockdown values ##
# need to remake at sgRNA level
features <- read.table('/raleighlab/data1/liuj/schwannoma/perturb-seq/scripts/libraries/features_perturb_1.csv',sep=',',header=TRUE,row.names=1)

# Normalize data to CPM
dat.sub.kd <- NormalizeData(dat.all, normalization.method = "RC", scale.factor = 1e6, verbose = TRUE) # CPM normalize

# build non targeting expression matrix 
targets <- sapply(strsplit(as.character(features$name),"_"), `[`, 1)
targets <- unique(targets)
targets <- targets[1:(length(targets)-2)]

sgRNATargetTbl <- sgRNATblConcordant[,c("feature_call","Gene")]
sgRNATargetTbl <- sgRNATargetTbl[!duplicated(sgRNATargetTbl$feature_call),]
row.names(sgRNATargetTbl) <- as.character(sgRNATargetTbl$feature_call)
sgRNATargetTbl <- sgRNATargetTbl[1:110,]

ntTargetTbl <- data.frame(row.names = row.names(sgRNATargetTbl), matrix(nrow = nrow(sgRNATargetTbl), ncol=length(smpid))) 
colnames(ntTargetTbl) <- smpid

# build bulk expression table for sgRNA on target cells 
kdTargetTbl <- data.frame(row.names = row.names(sgRNATargetTbl), matrix(nrow = nrow(sgRNATargetTbl), ncol=length(smpid))) 
colnames(kdTargetTbl) <- smpid

for (i in 1:length(smpid)){
  ntcells = sgRNATblConcordant[sgRNATblConcordant$.id == smpid[i],]
  ntcells = ntcells[complete.cases(ntcells),]
  ntcells = ntcells[grep("nc_",ntcells$feature_call),"cell_barcode"]
  
  if (length(ntcells) > 0 ){
    df = GetAssayData(dat.sub.kd, slot = 'data')
    colnames(df) = sapply(strsplit(colnames(df),"-"), `[`, 1)
    dfNt = df[intersect(targets,row.names(df)),intersect(ntcells,colnames(df))]
    ntTargetTbl[,i] <- apply(data.frame(dfNt),1,mean)[sgRNATargetTbl[row.names(ntTargetTbl),"Gene"]]
    
    for (j in 1:nrow(kdTargetTbl)){
      targetcells = sgRNATblConcordant[sgRNATblConcordant$.id == smpid[i] & sgRNATblConcordant$Gene == sgRNATargetTbl[row.names(kdTargetTbl[j,]),"Gene"],"cell_barcode"]
      if (length(targetcells) > 0 ){
        dfTarget = df[sgRNATargetTbl[row.names(kdTargetTbl[j,]),"Gene"],intersect(targetcells,colnames(df))]
        kdTargetTbl[j,i] <- mean(dfTarget,na.rm = TRUE)
      }
    }
  }
}

# get mean knockdown values - now use pseudocount 
remainingRNATbl <- round((kdTargetTbl+0.01)/(ntTargetTbl+0.01),3)

# generate ceiling of 1.0
remainingRNATbl[remainingRNATbl >= 1.0 ] = 1.0

# create heatmap of knockdown 
df <- melt(data.frame(gene=row.names(remainingRNATbl), remainingRNATbl))

plot <- ggplot(df, aes(x = gene,y=variable)) +
  geom_tile(aes(fill = value)) +
  geom_text(aes(label = round(value,2)),size=1) +
  scale_fill_gradient(low = "white", high = "steelblue",limit=c(0,1),na.value="gray50") +
  labs(title = "Mean sgRNA knockdown - pseudobulk", x="", y="") +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_blank(),panel.background=element_blank(),axis.line = element_line(size=0.25))
ggsave("pdx_perturb_meanKD_heatmap.png", plot, width =13, height = 1.8)

plot <- ggplot(df) +
  geom_histogram(aes(fill = variable,x=value), color=NA,binwidth =0.05,show.legend = FALSE) +
  scale_x_continuous(expand=c(0.01,0.01), limits = c(-0.05,1.05), breaks = seq(0,1,0.25)) +
  scale_y_continuous(expand=c(0.01,0.01), limits = c(0,70)) +
  labs(title = "Mean sgRNA knockdown - pseudobulk", x="mRNA remaining", y="Gene Targets-Condition Pair") +
  facet_wrap(~variable,scales="free", ncol = 4) +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=0,hjust=0.5,vjust=1),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_rect(colour='black',fill=NA, size=0.25),panel.background=element_blank(),axis.line = element_line(size=0.25))
ggsave("pdx_perturb_meanKD_histograms.png", plot, width =6, height = 2)


###################################################
# subset cells - get cell names with barcode suffixes and use them as row names
# USE ONLY single sgRNA CELLS 
smpidNumTbl <- data.frame(smpid = smpid, num = seq(1,length(smpid),1))
row.names(smpidNumTbl) <- as.character(smpidNumTbl$smpid)
sgRNATblConcordant$full_barcode <- paste(sgRNATblConcordant$cell_barcode,"-1_",smpidNumTbl[as.character(sgRNATblConcordant$.id),"num"],sep='')
row.names(sgRNATblConcordant) <- as.character(sgRNATblConcordant$full_barcode)

# subset to only wells with one sgRNA
dat.sub.vsinvitro <- subset(dat.sub.kd,cells = row.names(sgRNATblConcordant))
dat.sub.vsinvitro$sgRNA <- as.character(sgRNATblConcordant[colnames(dat.sub.vsinvitro),"feature_call"])
dat.sub.vsinvitro$sgRNA <- paste(dat.sub.vsinvitro$sgRNA, dat.sub.vsinvitro$RT,sep='_')

# get hit classes for metadata
sgRNAClass <- features[,c("name","target_gene_name")]
sgRNAClassRT <- rbind(data.frame(row.names = paste(row.names(sgRNAClass),"_0Gy",sep=""), sgRNAClass, cond=rep("0Gy",nrow(sgRNAClass))),
                      data.frame(row.names = paste(row.names(sgRNAClass),"_9Gy",sep=""), sgRNAClass, cond=rep("9Gy",nrow(sgRNAClass))),
                      data.frame(row.names = paste(row.names(sgRNAClass),"_12_5Gy",sep=""), sgRNAClass, cond=rep("12_5Gy",nrow(sgRNAClass))))
sgRNAClassRT <- sgRNAClassRT[,-1]

# GEM group normalize to non-targetings then implement David Wu's density UMAPs on only sgRNA expressing cells. Then do the modules
dfT.list <- list()

for (i in 1:length(smpid)){
  ntcells = sgRNATblConcordant[sgRNATblConcordant$.id == smpid[i],] # isolate each GEM group 
  ntcells = ntcells[complete.cases(ntcells),]
  tcells = ntcells[grep("nc_",ntcells$feature_call,invert=TRUE),"full_barcode"] # barcodes for on target cells 
  ntcells = ntcells[grep("nc_",ntcells$feature_call),"full_barcode"] # barcodes for non targeting cells 
  
  # get NT normalization vector 
  if (length(ntcells) > 0 ){
    df = GetAssayData(subset(dat.sub.vsinvitro, idents = smpid[i]), slot = 'data')
    #colnames(df) = sapply(strsplit(colnames(df),"-"), `[`, 1)
    dfNt = df[,intersect(ntcells,colnames(df))]
    ntVector<- apply(data.frame(dfNt),1,mean)
    
    dfT = round((df+0.01)/(ntVector+0.01),5)
    dfT.list[[i]] <- dfT
  }
}

names(dfT.list) <- smpid

dfT.gemnorm <- do.call(cbind.data.frame, dfT.list)
colnames(dfT.gemnorm) <- gsub("VS_perturb_1.","",colnames(dfT.gemnorm))
colnames(dfT.gemnorm) <- gsub("VS_perturb_2.","",colnames(dfT.gemnorm))
colnames(dfT.gemnorm) <- gsub("VS_perturb_3.","",colnames(dfT.gemnorm))
colnames(dfT.gemnorm) <- gsub("VS_perturb_4.","",colnames(dfT.gemnorm))

##############################################
# Normalize to Average of 0Gy both reps only #
##############################################

dfT0GyNorm.list <- list()

# generate pseudobulk NT vector for all
ntcells = sgRNATblConcordant[sgRNATblConcordant$cond == "0Gy",] # isolate each GEM group 
ntcells = ntcells[complete.cases(ntcells),]
ntcells = ntcells[grep("nc_",ntcells$feature_call),"full_barcode"]
df = GetAssayData(dat.sub.vsinvitro, slot = 'data')
dfNt = df[,intersect(ntcells,colnames(df))]
ntVector<- apply(data.frame(dfNt),1,mean)

for (i in 1:length(smpid)){
  df = GetAssayData(subset(dat.sub.vsinvitro, idents = smpid[i]), slot = 'data')
  dfT = round((df+0.01)/(ntVector+0.01),5)
  dfT0GyNorm.list[[i]] <- dfT
}

names(dfT0GyNorm.list) <- smpid

dfT0GyNorm <- do.call(cbind.data.frame, dfT0GyNorm.list)
colnames(dfT0GyNorm) <- gsub("VS_perturb_1.","",colnames(dfT0GyNorm))
colnames(dfT0GyNorm) <- gsub("VS_perturb_2.","",colnames(dfT0GyNorm))
colnames(dfT0GyNorm) <- gsub("VS_perturb_3.","",colnames(dfT0GyNorm))
colnames(dfT0GyNorm) <- gsub("VS_perturb_4.","",colnames(dfT0GyNorm))

# Identify only perturbations with KD
sgRNAWithKD <- remainingRNATbl[apply(remainingRNATbl[,1:2],1,min) < 0.25,] # at least one of the no RT reps have 
sgRNAWithKD <- row.names(sgRNAWithKD[!is.na(sgRNAWithKD$VS_perturb_1) & !is.na(sgRNAWithKD$VS_perturb_2),])
sgRNAWithKD <- c(sgRNAWithKD,"nc_01688","nc_01736","nc_02727","nc_03150")

# Expand sgRNA level metadata 
sgRNATblConcordant$cond <- sgRNATblConcordant$.id 
sgRNATblConcordant$cond <- gsub("VS_perturb_1","0Gy",sgRNATblConcordant$cond)
sgRNATblConcordant$cond <- gsub("VS_perturb_2","0Gy",sgRNATblConcordant$cond)
sgRNATblConcordant$cond <- gsub("VS_perturb_3","9Gy",sgRNATblConcordant$cond)
sgRNATblConcordant$cond <- gsub("VS_perturb_4","12_5Gy",sgRNATblConcordant$cond)
sgRNATblConcordant$geneCond <- paste(sgRNATblConcordant$Gene, sgRNATblConcordant$cond, sep=":")
sgRNATblConcordant$sgRNAcond <- paste(sgRNATblConcordant$feature_call, sgRNATblConcordant$cond, sep=":")

# export this table for use
write.table(sgRNATblConcordant,'vs_perturb_reanalysis_102022_sgRNATblConcordant.txt',sep='\t',quote=FALSE,row.names=TRUE,col.names=NA)

#######################################################################
# get scores for signatures from single cell HEI-193 RNA-seq clusters #
top50 <- read.table('/raleighlab/data1/liuj/schwannoma/scrna_cells/analysis/plots/top50markers.txt',header=TRUE,sep='\t') # Tim's scRNA data only
#top50 <- read.table('/raleighlab/data1/liuj/schwannoma/scrna_cells/analysis/combined/vs_perturb_scRNA_cells_top50markers_harmony_integration_res0.4.txt',header=TRUE,sep='\t') # integrated with perturb-seq
load('/raleighlab/data1/liuj/schwannoma/scrna/analysis_scrna/top50markers_20221107.Rds') # human schwannoma cell states

coi <- sgRNATblConcordant[matchAll(sgRNAWithKD,sgRNATblConcordant$feature_call),]

# use sgRNA level or gene level with radiation condition
sigRes <- data.frame(matrix(data = 0,nrow = length(unique(coi$geneCond)), ncol = 15))
row.names(sigRes) <- unique(coi$geneCond)
colnames(sigRes) <- paste("Signature",seq(0,14,1))
# colnames(sigRes) <- c("Immunogenic schwannoma", "Myelin remodeling schwannoma", "Microglia", "Vascular remodeling schwannoma ", "T cells", "Neural regeneration schwannoma",
#                       "Axon injury schwannoma", "Vascular endothelia", "Macrophages", "Complement activation schwannoma", "Stromal cells", "Progenitor schwannoma", "Vascular smooth muscle")
colnames(sigRes) <- c("Fe Metabolism","S phase","Gluthathione","EMT","Cell Stress (ATF3, TK1)","G2/M","TGFb","Cell Stress (ATF3, SGK1)","RNA splicing","Mitochondria",
                       "Ribosome","Cell membrane","p53 pathway","Interferon","Ribosome PLP2+")
#colnames(sigRes) <- c("S Phase", "Unfolded Protein Response", "G2-M Phase", "Hypoxia", "Pyrimidine Metabolism", 
#                      "p53 Pathway", "Integrin Signaling", "Ribosome 1", "Ribosome 2", "Cell membrane", "Translation factors",
#                      "Apoptosis", "Glycolysis", "RNA Splicing", "Interferon Signaling")
 
for (i in 1:nrow(sigRes)){
  for (j in 1:ncol(sigRes)){
    #goi = pull(top50[top50$cluster == j-1,"gene"][1:10,]) # top n marker genes for clinical samples
    goi = top50[top50$cluster == j-1,"gene"][1:10]
    #goi = goi[grep("AC092069.1",invert = TRUE,goi)]
    cells = row.names(sigRes)[i]
    
    # establish cutoff of needing at least x cells 
    if (nrow(sgRNATblConcordant[sgRNATblConcordant$geneCond == cells,]) >= 10) {
      #pseudobulk = apply(dfT.gemnorm[,row.names(sgRNATblConcordant[sgRNATblConcordant$geneCond == cells,])],1,mean) # pseudobulk using mean across all sgRNA conds 
      pseudobulk = apply(dfT0GyNorm[,row.names(sgRNATblConcordant[sgRNATblConcordant$geneCond == cells,])],1,mean) # pseudobulk using mean across all sgRNA conds, using 0 Gy normalized NTs 
      score = mean(pseudobulk[goi],na.rm = TRUE) # signature using mean 
    } else {
      score = NA
    }
    
    sigRes[i,j] <- score
  }
}

sigRes <- sigRes[complete.cases(sigRes),]
sigRes <- sigRes[order(row.names(sigRes)),]

#sigResClinical <- sigRes[,grep("schwannoma",colnames(sigRes))] # pull out schwannoma cell states only from clinical data 

sigRes <- cbind(sigRes,sigResClinical)
  
sgRNAClassRT <- data.frame(row.names = row.names(sigRes), genecond = row.names(sigRes), cond = str_split_fixed(row.names(sigRes), ":", 2)[,2])
statesClass <- data.frame(row.names = colnames(sigRes), cluster = colnames(sigRes),type = c(rep("states",15),rep("types",7)))


#####################
# cluster by groups #
#pdf("vs_perturb_modules_intHarmony_res04_GemNorm_split.pdf",width = 14, height=6)
#pdf("vs_perturb_modules_intHarmony_res04_0GyNorm_split.pdf",width = 14, height=6)
#pdf("vs_perturb_modules_res03_0GyNorm_split.pdf",width = 14, height=6)
#pdf("vs_perturb_modules_res03_GemNorm_split.pdf",width = 14, height=6)
#pdf("vs_perturb_clinical_modules_0GyNorm_split.pdf",width = 14, height=6)
pdf("vs_perturb_modules_res03_and_clinical_0GyNorm_split.pdf",width = 14, height=8)


Heatmap(as.matrix(t(log(sigRes[complete.cases(sigRes),],2))),
        top_annotation = HeatmapAnnotation(cond = sgRNAClassRT$cond,col = list(cond = c("0Gy" = "#1B9E77", "9Gy" = "#D95F02", "12_5Gy" = "#7570B3"))),
        column_split = factor(sgRNAClassRT$cond, levels = c("0Gy","9Gy","12_5Gy")),
        row_split = statesClass$type,
        cluster_column_slice = FALSE,
        width = ncol(t(sigRes))*unit(4, "mm"),
        height = nrow(t(sigRes))*unit(4, "mm"),
        col = colorRamp2(seq(-3, 3, length = 3), c("#2268AD", "#EEEEEE", "#B31B2C"))
        )
dev.off()



# deconvolute cell death signatures - AC092069.1 was driving a huge part of the ribosome and cell membrane clusters
goi = top50[top50$cluster == 14,"gene"][1:10]
pheatmap(dfT.gemnorm[goi,row.names(sgRNATblConcordant[sgRNATblConcordant$geneCond == "PTPRG:12_5Gy",])], cluster_rows = FALSE, cluster_cols = FALSE)

# generate heatmap using all marker genes.
df <- t(dfT0GyNorm[as.character(unique(top50[grep(),gene])),])
#df <- t(dfT.gemnorm[as.character(unique(top50$gene)),])
df <- data.frame(geneCond = sgRNATblConcordant[row.names(df),"geneCond"], df)
df <- aggregate(. ~ geneCond, data = df, FUN = mean, na.rm = TRUE)
row.names(df) <- df$geneCond
df <- df[,2:ncol(df)]
df <- t(df)
df <- df[,row.names(sigRes)]


pdf("vs_perturb_markerGenes_intHarmony_res04_0GyNorm.pdf",width = 16, height=80)
pheatmap(log(df[complete.cases(df),],2),border_color = NA,#annotation_col = sgRNAClassRT[,c("cond","cells")],
         cellwidth = 10,cellheight = 10,breaks = seq(-1,1,0.04),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color=rev(brewer.rdbu(50)))
dev.off()

#####################################
# Calculate differential expression #
dat.sub.vsinvitro$sgRNA <- as.character(sgRNATblConcordant[colnames(dat.sub.vsinvitro),"feature_call"])

# 0 Gy samples first
degResinvitro.list<-list()
allNTs = sgRNATblConcordant[sgRNATblConcordant$geneCond == "nc:0Gy",]

for (i in 1:length(sgRNAWithKD)){
  cells = dat.sub.vsinvitro
  cells = subset(cells, idents=c("VS_perturb_1","VS_perturb_2"))
  Idents(cells) <- 'sgRNA'
  cells =  subset(cells, idents = c(sgRNAWithKD[i],allNTs$feature_call))
  cells[["RNA"]]@counts <- as.matrix(cells[["RNA"]]@counts)+1 # add pseudocount
  dnsmp = min(c(ncol(subset(cells, idents=sgRNAWithKD[i])), ncol(subset(cells, idents=allNTs$feature_call)))) # downsample cells to lesser of the two
  
  if(dnsmp > 10){
    res = FindMarkers(cells, ident.1 = sgRNAWithKD[i], ident.2 = allNTs$feature_call[grep(sgRNAWithKD[i],allNTs$feature_call,invert=TRUE)],test.use = 'MAST',assay = 'RNA', max.cells.per.ident = dnsmp)
    res$avg_log2FC = log((res$pct.1+0.001)/(res$pct.2+0.001),2) # manual log2fc 
    res = res[order(-res$avg_log2FC),]
    degResinvitro.list[[i]] <- res
  }
  else{
    degResinvitro.list[[i]] <- NA
  }
}

names(degResinvitro.list) <- paste(sgRNAWithKD,":0Gy",sep='')

# 9 Gy
degResinvitro9Gy.list<-list()
allNTs = sgRNATblConcordant[sgRNATblConcordant$geneCond == "nc:9Gy",]

for (i in 1:length(sgRNAWithKD)){
  cells = dat.sub.vsinvitro
  cells = subset(cells, idents=c("VS_perturb_3"))
  Idents(cells) <- 'sgRNA'
  cells =  subset(cells, idents = c(sgRNAWithKD[i],allNTs$feature_call))
  cells[["RNA"]]@counts <- as.matrix(cells[["RNA"]]@counts)+1 # add pseudocount
  dnsmp = min(c(ncol(subset(cells, idents=sgRNAWithKD[i])), ncol(subset(cells, idents=allNTs$feature_call)))) # downsample cells to lesser of the two
  
  if(dnsmp > 10){
    res = FindMarkers(cells, ident.1 = sgRNAWithKD[i], ident.2 = allNTs$feature_call[grep(sgRNAWithKD[i],allNTs$feature_call,invert=TRUE)],test.use = 'MAST',assay = 'RNA', max.cells.per.ident = dnsmp)
    res$avg_log2FC = log((res$pct.1+0.001)/(res$pct.2+0.001),2) # manual log2fc 
    res = res[order(-res$avg_log2FC),]
    degResinvitro9Gy.list[[i]] <- res
  }
  else{
    degResinvitro9Gy.list[[i]] <- NA
  }
}

names(degResinvitro9Gy.list) <- paste(sgRNAWithKD,":9Gy",sep='')

# 12.5 Gy
degResinvitro12_5Gy.list<-list()
allNTs = sgRNATblConcordant[sgRNATblConcordant$geneCond == "nc:12_5Gy",]

for (i in 1:length(sgRNAWithKD)){
  cells = dat.sub.vsinvitro
  cells = subset(cells, idents=c("VS_perturb_4"))
  Idents(cells) <- 'sgRNA'
  cells =  subset(cells, idents = c(sgRNAWithKD[i],allNTs$feature_call))
  cells[["RNA"]]@counts <- as.matrix(cells[["RNA"]]@counts)+1 # add pseudocount
  dnsmp = min(c(ncol(subset(cells, idents=sgRNAWithKD[i])), ncol(subset(cells, idents=allNTs$feature_call)))) # downsample cells to lesser of the two
  
  if(dnsmp > 10){
    res = FindMarkers(cells, ident.1 = sgRNAWithKD[i], ident.2 = allNTs$feature_call[grep(sgRNAWithKD[i],allNTs$feature_call,invert=TRUE)],test.use = 'MAST',assay = 'RNA', max.cells.per.ident = dnsmp)
    res$avg_log2FC = log((res$pct.1+0.001)/(res$pct.2+0.001),2) # manual log2fc 
    res = res[order(-res$avg_log2FC),]
    degResinvitro12_5Gy.list[[i]] <- res
  }
  else{
    degResinvitro12_5Gy.list[[i]] <- NA
 }
}

names(degResinvitro12_5Gy.list) <- paste(sgRNAWithKD,":12_5Gy",sep='')

# Combine DEG results and eliminate empty results 
degRes <- c(degResinvitro.list,degResinvitro9Gy.list,degResinvitro12_5Gy.list)
v <- as.logical(as.logical(lapply(degRes,function(x) nrow(x) >0)))
v[is.na(v)] <- FALSE
degRes <- degRes[v]

# get summary of KDs
degResSig <- list()
for (i in 1:length(degRes)){
  degResSig[[i]] <- degRes[[i]][degRes[[i]]$p_val < 0.01 & abs(degRes[[i]]$avg_log2FC) > 1,]
}

names(degResSig) <- names(degRes)
degResTbl <- data.frame(genecond = names(degResSig), DEGS = as.numeric(lapply(degResSig,nrow)))

degResTblCmp <- data.frame(row.names= degResTbl$genecond, sgRNA = str_split_fixed(degResTbl$genecond, ":",2)[,1], cond = str_split_fixed(degResTbl$genecond, ":",2)[,2], DEGS = degResTbl$DEGS)
degResTblCmp <- reshape(degResTblCmp, idvar = "sgRNA", timevar = "cond", direction = "wide")
row.names(degResTblCmp) <- as.character(degResTblCmp$sgRNA)
degResTblCmp$radEnrich <- degResTblCmp$DEGS.9Gy - degResTblCmp$DEGS.0Gy
degResTblCmp$radEnrich[is.na(degResTblCmp$radEnrich)] <- 0

# plot 0 Gy vs 9 Gy DEGs 
ggplot() +
  geom_point(data=degResTblCmp[degResTblCmp$radEnrich < 40,],aes(x = DEGS.0Gy, y = DEGS.9Gy), fill = "#7570B3", size=2, shape = 21) +
  geom_point(data=degResTblCmp[degResTblCmp$radEnrich > 40,],aes(x = DEGS.0Gy, y = DEGS.9Gy), fill = "#D95F02", size=2, shape = 21) +
  geom_abline(slope = 1,intercept = 0, linetype='dashed') +
  geom_text_repel(data= degResTblCmp[degResTblCmp$radEnrich > 40,], aes(x = DEGS.0Gy, y = DEGS.9Gy, label=sgRNA),xnudge=0,ynudge=0,max.overlaps=25,size=2) +
  scale_x_continuous(limits = c(0,250)) +
  scale_y_continuous(limits = c(0,250)) +
  labs(title = "Number of DEGs using MAST\n(p < 0.05 and |log2FC| > 1)", x="DEG's 0 Gy", y="DEG's 5x 1.8 Gy") +
  theme_base() +
  theme(axis.text=element_text(size=10),text=element_text(size=10),axis.text.x=element_text(angle=0,hjust=0,vjust=1),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_blank(),panel.background=element_blank(),axis.line = element_line(size=0.25))
ggsave("vs_perturb_DEGs_scatter_MAST.pdf", width =5, height = 5) 

# export points 
write.table(degResTblCmp,'vs_perturb_reanalysis_102022_DEGs_scatter.txt',sep='\t',quote=FALSE,row.names=TRUE,col.names=NA)

#####################################################################
# DEGs and volcano plot for PTPRG, SOX6, APOC1, SHH? At gene levels #
cells = dat.sub.vsinvitro
cells$target = as.character(sgRNATblConcordant[colnames(cells),"Gene"])
cells$geneCond = as.character(sgRNATblConcordant[colnames(cells),"geneCond"])
cells[["RNA"]]@counts <- as.matrix(cells[["RNA"]]@counts)+1 # add pseudocount
Idents(cells) <- "geneCond"

goi <- "PTPRG"
dose <- ":9Gy"

# DEG in noRT conditions
dnsmp = min(c(ncol(subset(cells, idents=paste(goi,dose,sep=''))), ncol(subset(cells, idents=paste("nc",dose,sep=''))))) # downsample cells to lesser of the two
res = FindMarkers(cells, ident.1 = paste(goi,dose,sep=''), ident.2 = paste("nc",dose,sep=''),test.use = 'MAST',assay = 'RNA', max.cells.per.ident = dnsmp)
res$avg_log2FC = log((res$pct.1+0.001)/(res$pct.2+0.001),2) # manual log2fc 
res = res[order(-res$avg_log2FC),]
res$gene = row.names(res)
write.table(res,paste("vs_DE_gene_",goi,dose,".txt",sep=''),row.names=TRUE,col.names=NA,sep='\t',quote=FALSE)

ggplot() +
  geom_point(data=res[res$p_val >= 0.05 | abs(res$avg_log2FC) < 1,],aes(x = avg_log2FC, y= -log(p_val,10)), color = "gray50", size=2, shape = 19) +
  geom_point(data=res[res$p_val < 0.05 & res$avg_log2FC > 1,],aes(x = avg_log2FC, y= -log(p_val,10)), fill = "red", size=2, shape = 21) +
  geom_point(data=res[res$p_val < 0.05 & res$avg_log2FC < -1,],aes(x = avg_log2FC, y= -log(p_val,10)), fill = "blue", size=2, shape = 21) +
  geom_text_repel(data= head(res[res$p_val < 0.05 & res$avg_log2FC > 1,],15), aes(x = avg_log2FC, y = -log(p_val,10), label=gene),max.overlaps=25,size=2,segment.size=0.25) +
  geom_text_repel(data= tail(res[res$p_val < 0.05 & res$avg_log2FC < -1,],15), aes(x = avg_log2FC, y = -log(p_val,10), label=gene),max.overlaps=50,size=2,segment.size=0.25) +
  scale_x_continuous(limits = c(-8,8),breaks = seq(-8,8,1)) +
  scale_y_continuous(limits = c(0,7),breaks = seq(0,7,1)) +
  labs(title = paste("DEGs MAST (p < 0.05 and |log2FC| > 1)\n",goi,dose,sep=''), x="log2(fold change)\nsgRNA vs sgNT", y="-log10(pval)") +
  theme_base() +
  theme(axis.text=element_text(size=10),text=element_text(size=10),axis.text.x=element_text(angle=0,hjust=0.5,vjust=1),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_blank(),panel.background=element_blank(),axis.line = element_line(size=0.25))
ggsave(paste("vs_DE_gene_",goi,dose,"_scatter.pdf",sep=''), width =5, height = 5) 


############################################
# QC plots, sgRNA coverage, and knockdown  #
dat.sub.vsinvitro$sgRNA <- as.character(sgRNATblConcordant[colnames(dat.sub.vsinvitro),"feature_call"])
dat.sub.vsinvitro$target <- as.character(sgRNATblConcordant[colnames(dat.sub.vsinvitro),"Gene"])
Idents(dat.sub.vsinvitro) <- dat.sub.vsinvitro$sgRNA
dat.sub.mix <- subset(dat.sub.vsinvitro, idents = sgRNAWithKD) 

Idents(dat.sub.mix) <- "orig.ident"
dat.sub.mix$orig.ident <- gsub("VS_perturb_1","0Gy_1",dat.sub.mix$orig.ident)
dat.sub.mix$orig.ident <- gsub("VS_perturb_2","0Gy_2",dat.sub.mix$orig.ident)
dat.sub.mix$orig.ident <- gsub("VS_perturb_3","1.8Gyx5",dat.sub.mix$orig.ident)
dat.sub.mix$orig.ident <- gsub("VS_perturb_4","12.5Gyx1",dat.sub.mix$orig.ident)

SCpubr::do_ViolinPlot(sample = dat.sub.mix, features = c("nCount_RNA", "nFeature_RNA","percent.mt"),line_width =0.25,boxplot_width = 0.1,font.size = 6)
ggsave(filename = "vs_perturb_keep_QC_qc_violin.pdf",width = 4, height=2)

df <- data.frame(dat.sub.mix$orig.ident,dat.sub.mix$nCount_RNA,dat.sub.mix$nFeature_RNA,dat.sub.mix$percent.mt)
write.table(df,"vs_perturb_keep_QC_qc_violin.txt",quote=FALSE,sep='\t')

# get summary data - only for those with KD. 
df <- sgRNATblConcordant[colnames(dat.sub.mix),]
df <- melt(table(df[,c("Gene","cond","feature_call")]))

plot <- ggplot(df) +
  geom_bar(stat = 'identity',aes(x=Gene,y=value,fill=cond),position='dodge') +
  labs(title = "Number of cells with sgRNAs detected", x="", y="") +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_blank(),panel.background=element_blank(),axis.line = element_line(size=0.25))
ggsave("vs_perturb_keep_sgRNAs_hist.pdf", plot, width =3.5, height = 2)

write.table(df,"vs_perturb_keep_sgRNAs_hist.txt",quote=FALSE,sep='\t')

# KD histogram
# create heatmap of knockdown 
df <- data.frame(sgRNA = row.names(remainingRNATbl), remainingRNATbl)
df <- melt(df[as.character(unique(sgRNATblConcordant[colnames(dat.sub.mix),"feature_call"])),])
df <- df[!is.na(df$sgRNA),]

plot <- ggplot(df) +
  geom_histogram(aes(fill = variable,x=value), color=NA,binwidth =0.05,show.legend = FALSE) +
  scale_x_continuous(expand=c(0.01,0.01), limits = c(-0.05,1.05), breaks = seq(0,1,0.25)) +
  scale_y_continuous(expand=c(0.01,0.01), limits = c(0,25)) +
  labs(title = "Mean sgRNA knockdown - pseudobulk", x="mRNA remaining", y="Gene Targets-Condition Pair") +
  facet_wrap(~variable,scales="free", ncol = 4) +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=0,hjust=0.5,vjust=1),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_rect(colour='black',fill=NA, size=0.25),panel.background=element_blank(),axis.line = element_line(size=0.25))
ggsave("vs_perturb_keep_meanKD_histograms.pdf", plot, width = 6, height = 2)
write.table(df,"vs_perturb_keep_meanKD_histograms.txt",quote=FALSE,sep='\t')

# minimum KD histogram
df <- data.frame(sgRNA = row.names(remainingRNATbl), remainingRNATbl)
df <- df[as.character(unique(sgRNATblConcordant[colnames(dat.sub.mix),"feature_call"])),]
df <- data.frame(sgRNA = df$sgRNA, value = apply(df[,2:5],1,min))
df <- df[!is.na(df$sgRNA),]

plot <- ggplot(df) +
  geom_bar(stat = 'identity',aes(x=sgRNA,y=value), color=NA) +
  scale_y_continuous(expand=c(0.01,0.01), limits = c(-0.01,1.01), breaks = seq(0,1,0.25)) +
  labs(title = "Best sgRNA knockdown - pseudobulk", y="mRNA remaining", y="sgRNA") +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_rect(colour='black',fill=NA, size=0.25),panel.background=element_blank(),axis.line = element_line(size=0.25))
ggsave("vs_perturb_keep_meanKD_sgRNA_min.pdf", plot, width = 4, height = 2)
write.table(df,"vs_perturb_keep_meanKD_sgRNA_min.txt",quote=FALSE,sep='\t')

# run umap on subset of kept cells 
dat.sub.mix <- SCTransform(dat.sub.mix, verbose = TRUE)
dat.sub.mix <- RunPCA(dat.sub.mix, verbose = TRUE)

# PCA dimensions
ElbowPlot(object = dat.sub.mix,ndims = 50) 

# Run UMAP embeddings
mindist=0.5
dat.sub.mix <- RunUMAP(dat.sub.mix, reduction = "pca",dims = 1:30, min.dist = mindist)
DimPlot(object = dat.sub.mix, group.by = 'orig.ident', reduction = "umap", pt.size=0.3,raster=FALSE)
ggsave(filename = "vs_perturb_keep_umap_origidents.pdf",width = 4, height=3)

# cell cycle scoring 
dat.sub.mix <- CellCycleScoring(dat.sub.mix, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
SCpubr::do_FeaturePlot(dat.sub.mix, features = c("S.Score","G2M.Score"),order = TRUE)
ggsave(filename = "vs_perturb_keep_cellcycle.pdf",width = 7, height=5)

# embeddings for further exploration
embeddings <- Embeddings(dat.sub.mix, reduction="umap")
embeddings <- data.frame(embeddings, orig.ident=as.character((dat.sub.mix$orig.ident[row.names(embeddings)])),
                         S.Score=as.character((dat.sub.mix$S.Score[row.names(embeddings)])),
                         G2M.Score=as.character((dat.sub.mix$G2M.Score[row.names(embeddings)])))
write.table(embeddings,'vs_perturb_keep_umap_origidents_cc.txt',sep='\t',quote=FALSE,row.names=TRUE,col.names=NA)

# output passed sgRNA metadata 
features.sub <- features[names(table(dat.sub.mix$sgRNA)),]
write.table(features.sub,'vs_perturb_keep_featuresSub.txt',sep='\t',quote=FALSE,row.names=TRUE,col.names=NA)
