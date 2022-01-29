##Seurat analysis timer series dataset hiPSC-SANCM
##this script was run using Seurat v3.2.1 
#remotes::install_version("Seurat", version = "3.2.1")
#the following packages were used

library(Seurat)
library(URD)
library(Matrix)
library(fossil) 
library(dplyr)
library(plyr)
library(ggplot2)
library(cowplot)
library(patchwork)


#load seurat object
setwd("/path/to/directory/")
data<-readRDS("210705_hiPSC_SANCM.rds")

data
#An object of class Seurat 
#18837 features across 4087 samples within 1 assay 
#Active assay: RNA (18837 features, 0 variable features)

##QC
#calculate percent mitochondrial reads
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#and as scatter plots
plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1
plot2
#filter dataset to remove dying cells and potential doublets
data <- subset(data, subset = nFeature_RNA > 600 & nCount_RNA < 100000 & percent.mt < 50)
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

data
#An object of class Seurat 
#18834 features across 3300 samples within 1 assay 
#Active assay: RNA (18834 features, 0 variable features)


##SCTranform normalization + batch correction
#split the dataset into a list of 3 seurat objects (plate/batch 1-3)
data.list <- SplitObject(data, split.by = "batch")

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
#we set normalization.method = 'SCT'
data.anchors <- FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT", 
                                       anchor.features = data.features, verbose = TRUE)
data.integrated <- IntegrateData(anchorset = data.anchors, normalization.method = "SCT",verbose = TRUE)

#we continue working in the assay "integrated"
DefaultAssay(data.integrated) <- "integrated"

##dimensionality reduction
#proceed with downstream analysis (i.e. visualization, clustering) on the integrated dataset
data.integrated <- RunPCA(data.integrated, verbose = FALSE)
##Elbowplot was additionally used to determine dimensionality of dataset
ElbowPlot(data.integrated, ndims = 50)

#the top 20 PC were chosen to cluster cells
data.integrated <- FindNeighbors(data.integrated, dims = 1:20)
data.integrated <- FindClusters(object = data.integrated, reduction = "pca",dims = 1:20, resolution = 0.6,random.seed = 2020)
head(Idents(data.integrated), 5)
#UMAP (as shown in supplementary Fig3B)
data.integrated <- RunUMAP(data.integrated, dims = 1:25)
DimPlot(data.integrated, group.by = "orig.ident", label = TRUE, repel = TRUE)

##differential gene expression
#for differential gene expression after integration, we switch back to "RNA" assay and perform logNorm
DefaultAssay(data.integrated) <- "RNA"
# Normalize RNA data for visualization purposes
data.integrated <- NormalizeData(data.integrated, verbose = FALSE)
data.integrated <- ScaleData(data.integrated, verbose = FALSE)

#Heatmap
all.markers <- FindAllMarkers(object = data.integrated, only.pos = TRUE, min.pct = 0.25,  thresh.use = 0.25)
top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(object = data.integrated, features  = top10$gene, label = TRUE)

#Heatmap shows clusters enriched in ERCC-spike in DNA, which can also be identified as cells with very low counts
VlnPlot(data.integrated, "nFeature_RNA")

#cluster 7 and 13 were highly enriched in spike-in DNA and showed very low count of genes
#therefore they were removed for further analysis
subset <- SubsetData(object = data.integrated, ident.remove = c("7", "13"))

#cluster annotation according to the collected timepoint
current.cluster.ids <- c(0,1,2,3,4,5,6,8,9,10,11,12)
new.cluster.ids <- c("D10_2", "D06", "D19_1", "D05", "D0_1", "D04_1", "D0_2",
                     "D10_1", "D19_2", "D06", "D04_2", "D19_3")
subset@active.ident <- plyr::mapvalues(x = subset@active.ident, from = current.cluster.ids, to = new.cluster.ids)
#and ordered accordingly
my_levels <- c("D0_1", "D0_2", "D04_1", "D04_2", "D05", "D06", "D10_1", "D10_2", "D19_1", "D19_2", "D19_3")
subset@active.ident <- factor(x = subset@active.ident, levels = my_levels)
DimPlot(subset, reduction = "umap", label = TRUE, repel = TRUE)

#caluclate DE for each cluster (Supplementary TableS6)
#markers of D0_1
markers_D0_1 <- FindMarkers(subset, ident.1 = "D0_1", ident.2 = NULL, only.pos = TRUE) 
markers_D0_1[
  with(markers_D0_1, order( avg_logFC, decreasing = TRUE)),
]
write.csv(markers_D0_1, file="/path/to/directory/markers_D0_1.csv", row.names = TRUE)

#markers of D0_2
markers_D0_2 <- FindMarkers(subset, ident.1 = "D0_2", ident.2 = NULL, only.pos = TRUE) 
markers_D0_2[
  with(markers_D0_2, order( avg_logFC, decreasing = TRUE)),
]
write.csv(markers_D0_2, file="/path/to/directory/markers_D0_2.csv", row.names = TRUE)

#markers of D04_1
markers_D04_1 <- FindMarkers(subset, ident.1 = "D04_1", ident.2 = NULL, only.pos = TRUE) 
markers_D04_1[
  with(markers_D04_1, order( avg_logFC, decreasing = TRUE)),
]
write.csv(markers_D04_1, file="/path/to/directory/markers_D04_1.csv", row.names = TRUE)


#markers of D04_2
markers_D04_2 <- FindMarkers(subset, ident.1 = "D04_2", ident.2 = NULL, only.pos = TRUE) 
markers_D04_2[
  with(markers_D04_2, order( avg_logFC, decreasing = TRUE)),
]
write.csv(markers_D04_2, file="/path/to/directory/markers_D04_2.csv", row.names = TRUE)


#markers of D05
markers_D05 <- FindMarkers(subset, ident.1 = "D05", ident.2 = NULL, only.pos = TRUE) 
markers_D05[
  with(markers_D05, order( avg_logFC, decreasing = TRUE)),
]
write.csv(markers_D05, file="/path/to/directory/markers_D05.csv", row.names = TRUE)

#markers of D06
markers_D06 <- FindMarkers(subset, ident.1 = "D06", ident.2 = NULL, only.pos = TRUE) 
markers_D06[
  with(markers_D06, order( avg_logFC, decreasing = TRUE)),
]
write.csv(markers_D06, file="/path/to/directory/markers_D06.csv", row.names = TRUE)

#markers of D10_1
markers_D10_1 <- FindMarkers(subset, ident.1 = "D10_1", ident.2 = NULL, only.pos = TRUE) 
markers_D10_1[
  with(markers_D10_1, order( avg_logFC, decreasing = TRUE)),
]
write.csv(markers_D10_1, file="/path/to/directory/markers_D10_1.csv", row.names = TRUE)

#markers of D10_2
markers_D10_2 <- FindMarkers(subset, ident.1 = "D10_2", ident.2 = NULL, only.pos = TRUE) 
markers_D10_2[
  with(markers_D10_2, order( avg_logFC, decreasing = TRUE)),
]
write.csv(markers_D10_2, file="/path/to/directory/markers_D10_2.csv", row.names = TRUE)


#markers of D19_3
markers_D19_3 <- FindMarkers(subset, ident.1 = "D19_3", ident.2 = NULL, only.pos = TRUE) 
markers_D19_3[
  with(markers_D19_3, order( avg_logFC, decreasing = TRUE)),
]
write.csv(markers_D19_3, file="/path/to/directory/markers_D19_3.csv", row.names = TRUE)


#markers of D19_2
markers_D19_2 <- FindMarkers(subset, ident.1 = "D19_2", ident.2 = NULL, only.pos = TRUE) 
markers_D19_2[
  with(markers_D19_2, order( avg_logFC, decreasing = TRUE)),
]
write.csv(markers_D19_2, file="/path/to/directory/markers_D19_2.csv", row.names = TRUE)

#markers of D19_1
markers_D19_1 <- FindMarkers(subset, ident.1 = "D19_1", ident.2 = NULL, only.pos = TRUE) 
markers_D19_1[
  with(markers_D19_1, order( avg_logFC, decreasing = TRUE)),
]
write.csv(markers_D19_1, file="/path/to/directory/markers_D19_1.csv", row.names = TRUE)

##figures for paper
#UMAP Fig4B
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(8, "Spectral"))
colourCount = 11
DimPlot(subset, reduction = "umap", label = TRUE, pt.size = 0.8, cols = getPalette(colourCount))

#Heatmap Supplementary FigS3A
all.markers <- FindAllMarkers(object = subset, only.pos = TRUE, min.pct = 0.25,  thresh.use = 0.25)
top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(object = subset, features  = top10$gene, label = T, group.colors = getPalette(colourCount)) + scale_fill_gradientn(colors =rev(brewer.pal(n = 11, name = "RdBu")))

#FeaturePlots and ViolinPlot Fig4C
FeaturePlot(subset, features = c("TNNT2", "ACTN2"), cols = rev(brewer.pal(n = 11, name = "RdBu")))
VlnPlot(object = subset, features =c("TNNT2", "ACTN2"),cols = getPalette(colourCount))

#VlnPlots Fig4D-H
VlnPlot(object = subset, features =c("SOX2", "NANOG", "POU5F1", "EOMES", 
                                     "MESP1", "MESP2"),cols = getPalette(colourCount))
VlnPlot(object = subset, features =c("HOXA1", "NR2F2", "TBX5", "TBX18", "ISL1", "SHOX2", 
                                     "TBX3", "ALDH1A2", "PDPN", "WT1"),cols = getPalette(colourCount))
#ViolinPlot Supplementary FigS3C
VlnPlot(object = subset, features =c( "FOXA2", "SOX17"),cols = getPalette(colourCount))

##################################################################################################################################

#### SessionInfo ####
writeVersions <- function(sessionDir="/path/to/directory/Documentation"){
  write(paste0("Bioconductor version ", capture.output(tools:::.BioC_version_associated_with_R_version()),"\n"), 
        paste0(sessionDir,"/sessionInfo.txt"))
  write(capture.output(sessionInfo()), paste0(sessionDir,"/sessionInfo.txt"), append=TRUE)
}

#writeLines(capture.output(sessionInfo()), "../Documentation/sessionInfo.txt")
writeVersions()


##################################################################################################################################

