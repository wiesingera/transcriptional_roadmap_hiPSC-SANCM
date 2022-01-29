##Seurat analysis D19 endtimepoint hiPSC-SANCM and hiPSC-VCMs
##this script was run using Seurat v3.2.1 
#remotes::install_version("Seurat", version = "3.2.1")
#the following packages were loaded
library(Seurat)
library(Matrix)
library(fossil) 
library(dplyr)
library(plyr)
library(ggplot2)
library(cowplot)
library(patchwork)
library(RColorBrewer)

#load seurat object
setwd("/path/to/directory/")
D19<-readRDS("D19_SANCM_VCM_rep1_rep2.rds")

##QC
#calculate percent mitochondrial reads
D19[["percent.mt"]] <- PercentageFeatureSet(D19, pattern = "^MT-")
#Visualize QC metrics as a violin plot
VlnPlot(D19, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#and as scatter plots
plot1 <- FeatureScatter(D19, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(D19, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1+plot2
#filter dataset to remove dying cells and potential doublets
D19 <- subset(D19, subset = nFeature_RNA > 1000 & nFeature_RNA < 9000 & nCount_RNA < 60000 & percent.mt < 50)
VlnPlot(D19, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


##SCTranform normalization + batch correction
#split the dataset into a list of 3 seurat objects (plate/batch 1-3)
data.list <- SplitObject(D19, split.by = "batch")

#Prior to finding anchors, we perform standard preprocessing, normalization
#and identify variable features individually for each dataset. Variable feature selection is based on a 
#variance stabilizing transformation ("vst")
for (i in 1:length(data.list)) { data.list[[i]] <- SCTransform(data.list[[i]], verbose = TRUE)}

#select features for downstream integration, and run PrepSCTIntegration, which 
#ensures that all necessary Pearson residuals have been calculated
data.features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = 3000)
data.list <- PrepSCTIntegration(object.list = data.list, anchor.features = data.features, 
                                verbose = TRUE)

#identify anchors and integrate the datasets
#normalization.method = 'SCT'
data.anchors <- FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT", 
                                       anchor.features = data.features, verbose = TRUE)
data.integrated <- IntegrateData(anchorset = data.anchors, normalization.method = "SCT",verbose = TRUE)

#after integration, the default assay is automatically set to "integrated"
DefaultAssay(data.integrated) <- "integrated"

##dimensionality reduction
#proceed with downstream analysis (i.e. visualization, clustering) on the integrated dataset
#do not scale integrated data
data.integrated <- RunPCA(data.integrated, verbose = FALSE)
ElbowPlot(data.integrated, ndims = 50)

#cluster cells
data.integrated <- FindNeighbors(data.integrated, dims = 1:15)
data.integrated <- FindClusters(object = data.integrated, reduction = "pca",dims = 1:20, resolution = 0.6,random.seed = 2020)
head(Idents(data.integrated), 5)
#UMAP (as shown in supplementary Fig2A)
data.integrated <- RunUMAP(data.integrated, dims = 1:15)
DimPlot(data.integrated, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(data.integrated, group.by = "orig.ident", label = TRUE, repel = TRUE)

##differential gene expression
#for differential gene expression after integration, we switch back to "RNA" assay and perform logNorm
DefaultAssay(data.integrated) <- "RNA"
#Normalize original data for visualization purposes
data.integrated <- NormalizeData(data.integrated, verbose = FALSE)
data.integrated <- ScaleData(data.integrated, verbose = FALSE)

#Heatmap 
all.markers <- FindAllMarkers(object = data.integrated, only.pos = TRUE, min.pct = 0.25,  thresh.use = 0.25)
top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(object = data.integrated, features  = top10$gene, label = TRUE)

#Heatmap shows clustered enriched in ERCC-spike in DNA, which can also be identified as cells with very low counts
VlnPlot(data.integrated, "nFeature_RNA")
plot1 <- FeatureScatter(data.integrated, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(data.integrated, feature1 = "nFeature_RNA", feature2 = "nCount_RNA")
plot1
plot2

#clusters annotation: order of clusters were customized
current.cluster.ids <- c(0,1,2,3,4,5,6,7,8,9,10,11)
new.cluster.ids <- c("6", "1", "2", "5", "7", "4", "11", "3","9", "0", "8", "10")
data.integrated@active.ident <- plyr::mapvalues(x = data.integrated@active.ident, from = current.cluster.ids, to = new.cluster.ids)
#and ordered accordingly
my_levels <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11")
data.integrated@active.ident <- factor(x = data.integrated@active.ident, levels = my_levels)
DimPlot(data.integrated, reduction = "umap", label = TRUE, repel = TRUE)

#caluclate DE for each clusters (Supplementary TableS5)
#cluster0
markers_0 <- FindMarkers(data.integrated, ident.1 = "0", ident.2 = NULL, only.pos = TRUE) 
markers_0[
  with(markers_0, order( avg_logFC, decreasing = TRUE)),
]
write.csv(markers_0, file="/path/to/directory/markers_0.csv", row.names = TRUE)

#cluster1
markers_1 <- FindMarkers(data.integrated, ident.1 = "1", ident.2 = NULL, only.pos = TRUE) 
markers_1[
  with(markers1, order( avg_logFC, decreasing = TRUE)),
]
write.csv(markers_1, file="/path/to/directory/markers_1.csv", row.names = TRUE)

#cluster2
markers_2 <- FindMarkers(data.integrated, ident.1 = "2", ident.2 = NULL, only.pos = TRUE) 
markers_2[
  with(markers_2, order( avg_logFC, decreasing = TRUE)),
]
write.csv(markers_2, file="/path/to/directory/markers_2.csv", row.names = TRUE)

#cluster3
markers_3 <- FindMarkers(data.integrated, ident.1 = "3", ident.2 = NULL, only.pos = TRUE) 
markers_3[
  with(markers_3, order( avg_logFC, decreasing = TRUE)),
]
write.csv(markers_3, file="/path/to/directory/markers_3.csv", row.names = TRUE)

#cluster4
markers_4 <- FindMarkers(data.integrated, ident.1 = "4", ident.2 = NULL, only.pos = TRUE) 
markers_4[
  with(markers_4, order( avg_logFC, decreasing = TRUE)),
]
write.csv(markers_4, file="/path/to/directory/markers_4.csv", row.names = TRUE)

#cluster5
markers_5 <- FindMarkers(data.integrated, ident.1 = "5", ident.2 = NULL, only.pos = TRUE) 
markers_5[
  with(markers_5, order( avg_logFC, decreasing = TRUE)),
]
write.csv(markers_5, file="/path/to/directory/markers_5.csv", row.names = TRUE)

#cluster6
markers_6 <- FindMarkers(data.integrated, ident.1 = "6", ident.2 = NULL, only.pos = TRUE) 
markers_6[
  with(markers_6, order( avg_logFC, decreasing = TRUE)),
]
write.csv(markers_6, file="/path/to/directory/markers_6.csv", row.names = TRUE)

#cluster7
markers_7 <- FindMarkers(data.integrated, ident.1 = "7", ident.2 = NULL, only.pos = TRUE) 
markers_7[
  with(markers_7, order( avg_logFC, decreasing = TRUE)),
]
write.csv(markers_7, file="/path/to/directory/markers_7.csv", row.names = TRUE)

#cluster8
markers_8 <- FindMarkers(data.integrated, ident.1 = "8", ident.2 = NULL, only.pos = TRUE) 
markers_8[
  with(markers_8, order( avg_logFC, decreasing = TRUE)),
]
write.csv(markers_8, file="/path/to/directory/markers_8.csv", row.names = TRUE)

#cluster9
markers_9 <- FindMarkers(data.integrated, ident.1 = "9", ident.2 = NULL, only.pos = TRUE) 
markers_9[
  with(markers_9, order( avg_logFC, decreasing = TRUE)),
]
write.csv(markers_9, file="/path/to/directory/markers_9.csv", row.names = TRUE)

#cluster10
markers_10 <- FindMarkers(data.integrated, ident.1 = "10", ident.2 = NULL, only.pos = TRUE) 
markers_10[
  with(markers_10, order( avg_logFC, decreasing = TRUE)),
]
write.csv(markers_10, file="/path/to/directory/markers_10.csv", row.names = TRUE)

#cluster11
markers_11 <- FindMarkers(data.integrated, ident.1 = "11", ident.2 = NULL, only.pos = TRUE) 
markers_11[
  with(markers_11, order( avg_logFC, decreasing = TRUE)),
]
write.csv(markers_11, file="/path/to/directory/markers_11.csv", row.names = TRUE)


#cluster 9 are dead cells/empty wells as indicated by the enrichment in ERCC spike-in DNA, cluster 10 DE shows 
#exclusively cell cycle related genes and are excluded from further analysis. Additionally, cluster 11 could not be assigned
#to a specific cell type.Therefore, cluster 9, 10 and 11 were excluded from further analysis
subset <- SubsetData(object = data.integrated, ident.remove = c("9", "10", "11"))


##figures for paper
#colors for subset dataset
cols<-c("#993404","#fd8d3c", "#bd0026", "#fecc5c",  "#a1dab4", "#2b8cbe", "#016c59", "#bdd7e7", "#08519c")

#UMAP Fig3A
DimPlot(subset, reduction = "umap", label = TRUE, repel = TRUE, cols = cols)

#FeaturePlot and ViolinPlot Fig3B
FeaturePlot(subset, features = c("TNNT2", "ACTN2"), cols = rev(brewer.pal(n = 11, name = "RdBu")))
VlnPlot(object = subset, features =c("TNNT2", "ACTN2"), cols=cols)

#Heatmap Fig3C
all.markers <- FindAllMarkers(object = subset, only.pos = TRUE, min.pct = 0.25,  thresh.use = 0.25)
top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(object = subset, features  = top10$gene, label = T, group.colors = cols) + scale_fill_gradientn(colors =rev(brewer.pal(n = 11, name = "RdBu")))

#ViolinPlots Fig3D,E,F
VlnPlot(object = subset, features =c("TBX3", "ISL1", "BMP4", "SHOX2", "RGS6", "TBX18","NKX2-5", "CPNE5","WT1", "BNC1",), cols=cols)

##Supplementary Figures
#UMAP Supplementary Fig2A
DimPlot(data.integrated, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(data.integrated, group.by = "orig.ident", label = TRUE, repel = TRUE)

#Heatmap Supplementary Fig2B 
suppl <- SubsetData(data.integrated, idents = c("9", "10", "11"))
all.markers <- FindAllMarkers(object = suppl, only.pos = TRUE, min.pct = 0.25,  thresh.use = 0.25)
top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(object = suppl, features  = top10$gene, label = T, group.colors = cols) + scale_fill_gradientn(colors =rev(brewer.pal(n = 11, name = "RdBu")))

#ViolinPlots Supplementary Fig 2C-H
VlnPlot(object = subset, features =c("MYL2", "IRX4", "HOPX", "HEY2"), cols=cols)
VlnPlot(object = subset, features =c("HAPLN1", "GPC3", "SEMA3C", "PITX2"), cols=cols)
VlnPlot(object = subset, features =c("NFATC1", "FOXC1", "NRG1", "NPR3"), cols=cols)
VlnPlot(object = subset, features =c("MYH6", "HAMP", "NPPA"), cols=cols)
VlnPlot(object = subset, features =c("COL1A1","COL1A2", "COL3A1"), cols=cols)
VlnPlot(object = subset, features =c("VSNL1","GNAO1"), cols=cols)


#### SessionInfo ####
writeVersions <- function(sessionDir="/path/to/directory/Documentation"){
  write(paste0("Bioconductor version ", capture.output(tools:::.BioC_version_associated_with_R_version()),"\n"), 
        paste0(sessionDir,"/sessionInfo.txt"))
  write(capture.output(sessionInfo()), paste0(sessionDir,"/sessionInfo.txt"), append=TRUE)
}

#writeLines(capture.output(sessionInfo()), "../Documentation/sessionInfo.txt")
writeVersions()

