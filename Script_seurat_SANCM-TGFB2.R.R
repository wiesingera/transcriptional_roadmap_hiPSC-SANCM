library(Seurat)
library(Matrix)
library(fossil) 
library(dplyr)
library(plyr)
library(ggplot2)

setwd("/path/to/directory/")
data<-readRDS("data_TGFB2_SAN_ACM.rds")

data
#An object of class Seurat 
#18153 features across 1964 samples within 1 assay 
#Active assay: RNA (18153 features, 0 variable features)

## vizualize QC

VlnPlot(data, "nFeature_RNA")

# Visualize QC metrics as a violin plot
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1
plot2


#Normalize data

data.list <- SplitObject(data, split.by = "batch")

##Prior to finding anchors, we perform standard preprocessing, 
##and identify variable features individually for each. Note that Seurat v3 
##implements an improved method for variable feature selection based on a 
##variance stabilizing transformation ("vst")

for (i in 1:length(data.list)) { data.list[[i]] <- SCTransform(data.list[[i]], verbose = TRUE)}

#select features for downstream integration, and run PrepSCTIntegration, which 
#ensures that all necessary Pearson residuals have been calculated

data.features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = 3000)
data.list <- PrepSCTIntegration(object.list = data.list, anchor.features = data.features, 
                                verbose = TRUE)

#identify anchors and integrate the datasets
#make sure to set normalization.method = 'SCT'

data.anchors <- FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT", 
                                       anchor.features = data.features, verbose = TRUE)
data.integrated <- IntegrateData(anchorset = data.anchors, normalization.method = "SCT",verbose = TRUE)



##we continue working in the assay "integrated"
DefaultAssay(data.integrated) <- "integrated"

#proceed with downstream analysis (i.e. visualization, clustering) on the integrated dataset
#do not scale integrated data

data.integrated <- RunPCA(data.integrated, verbose = FALSE)
DimPlot(data.integrated, reduction = "pca", dims = c(1,2))
DimHeatmap(data.integrated, dims = 1, cells = 500, balanced = TRUE)

ElbowPlot(data.integrated, ndims = 50)

#15 dims were chosen
data.integrated <- RunUMAP(data.integrated, dims = 1:15)
DimPlot(data.integrated, reduction = "umap")
data.integrated <- FindNeighbors(data.integrated, dims = 1:15)
data.integrated <- FindClusters(object = data.integrated, reduction = "pca",dims = 1:15, resolution = 0.6,random.seed = 2020)
head(Idents(data.integrated), 5)
data.integrated <- RunUMAP(data.integrated, dims = 1:15)
DimPlot(data.integrated, reduction = "umap", label = TRUE, repel = TRUE)

DimPlot(data.integrated, group.by = "orig.ident", label = TRUE, repel = TRUE)

DefaultAssay(data.integrated) <- "RNA"


# Normalize RNA data for visualization purposes
data.integrated <- NormalizeData(data.integrated, verbose = FALSE)
data.integrated <- ScaleData(data.integrated, verbose = FALSE)

FeaturePlot(data.integrated, c("TNNT2", "SHOX2", "WT1", "NPPA", "KCNIP2", "HAMP", "ADM"))

all.markers <- FindAllMarkers(object = data.integrated, only.pos = TRUE, min.pct = 0.25,  thresh.use = 0.25)

top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

DoHeatmap(object = data.integrated, features  = top10$gene, label = TRUE)

VlnPlot(data.integrated, "nFeature_RNA")
plot1 <- FeatureScatter(data.integrated, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data.integrated, feature1 = "nFeature_RNA", feature2 = "nCount_RNA")
plot1
plot2

#remove cluster 6 (low-quality cells cells)
subset <- subset(data.integrated, idents = c("0","1", "2", "3", "4", "5", "7", "8", "9", "10"))


subset
#An object of class Seurat 
#38086 features across 1821 samples within 3 assays 
#Active assay: RNA (18153 features, 0 variable features)
# other assays present: SCT, integrated
#2 dimensional reductions calculated: pca, umap

DimPlot(subset, group.by = "orig.ident", label = TRUE, repel = TRUE) + NoLegend()
DimPlot(subset, reduction = "umap", label = TRUE, repel = TRUE)

#get maker genes of all clusters
##markers of 0
markers_0 <- FindMarkers(subset, ident.1 = "0", ident.2 = NULL, only.pos = TRUE) 
markers_0[
  with(markers_0, order( avg_log2FC, decreasing = TRUE)),
]

write.csv(markers_0, file="L:/basic/Personal Archive/A/awiesinger/scRNAseq_analyses/01Manuscript/TGFB2_scRNASeq/Analysis/DE_15dims/markers_0.csv", row.names = TRUE)

##markers of 1
markers_1 <- FindMarkers(subset, ident.1 = "1", ident.2 = NULL, only.pos = TRUE) 
markers_1[
  with(markers_1, order( avg_log2FC, decreasing = TRUE)),
]

write.csv(markers_1, file="L:/basic/Personal Archive/A/awiesinger/scRNAseq_analyses/01Manuscript/TGFB2_scRNASeq/Analysis/DE_15dims/markers_1.csv", row.names = TRUE)

##markers of 2
markers_2 <- FindMarkers(subset, ident.1 = "2", ident.2 = NULL, only.pos = TRUE) 
markers_2[
  with(markers_2, order( avg_log2FC, decreasing = TRUE)),
]

write.csv(markers_2, file="L:/basic/Personal Archive/A/awiesinger/scRNAseq_analyses/01Manuscript/TGFB2_scRNASeq/Analysis/DE_15dims/markers_2.csv", row.names = TRUE)

##markers of 3
markers_3 <- FindMarkers(subset, ident.1 = "3", ident.2 = NULL, only.pos = TRUE) 
markers_3[
  with(markers_3, order( avg_log2FC, decreasing = TRUE)),
]

write.csv(markers_3, file="L:/basic/Personal Archive/A/awiesinger/scRNAseq_analyses/01Manuscript/TGFB2_scRNASeq/Analysis/DE_15dims/markers_3.csv", row.names = TRUE)

##markers of 4
markers_4 <- FindMarkers(subset, ident.1 = "4", ident.2 = NULL, only.pos = TRUE) 
markers_4[
  with(markers_4, order( avg_log2FC, decreasing = TRUE)),
]

write.csv(markers_4, file="L:/basic/Personal Archive/A/awiesinger/scRNAseq_analyses/01Manuscript/TGFB2_scRNASeq/Analysis/DE_15dims/markers_4.csv", row.names = TRUE)

##markers of 5
markers_5 <- FindMarkers(subset, ident.1 = "5", ident.2 = NULL, only.pos = TRUE) 
markers_5[
  with(markers_5, order( avg_log2FC, decreasing = TRUE)),
]

write.csv(markers_5, file="L:/basic/Personal Archive/A/awiesinger/scRNAseq_analyses/01Manuscript/TGFB2_scRNASeq/Analysis/DE_15dims/markers_5.csv", row.names = TRUE)


##markers of 7
markers_7 <- FindMarkers(subset, ident.1 = "7", ident.2 = NULL, only.pos = TRUE) 
markers_7[
  with(markers_7, order( avg_log2FC, decreasing = TRUE)),
]

write.csv(markers_7, file="L:/basic/Personal Archive/A/awiesinger/scRNAseq_analyses/01Manuscript/TGFB2_scRNASeq/Analysis/DE_15dims/markers_7.csv", row.names = TRUE)

##markers of 8
markers_8 <- FindMarkers(subset, ident.1 = "8", ident.2 = NULL, only.pos = TRUE) 
markers_8[
  with(markers_8, order( avg_log2FC, decreasing = TRUE)),
]

write.csv(markers_8, file="L:/basic/Personal Archive/A/awiesinger/scRNAseq_analyses/01Manuscript/TGFB2_scRNASeq/Analysis/DE_15dims/markers_8.csv", row.names = TRUE)

##markers of 9
markers_9 <- FindMarkers(subset, ident.1 = "9", ident.2 = NULL, only.pos = TRUE) 
markers_9[
  with(markers_9, order( avg_log2FC, decreasing = TRUE)),
]

write.csv(markers_9, file="L:/basic/Personal Archive/A/awiesinger/scRNAseq_analyses/01Manuscript/TGFB2_scRNASeq/Analysis/DE_15dims/markers_9.csv", row.names = TRUE)

##markers of 10
markers_10 <- FindMarkers(subset, ident.1 = "10", ident.2 = NULL, only.pos = TRUE) 
markers_10[
  with(markers_10, order( avg_log2FC, decreasing = TRUE)),
]

write.csv(markers_10, file="L:/basic/Personal Archive/A/awiesinger/scRNAseq_analyses/01Manuscript/TGFB2_scRNASeq/Analysis/DE_15dims/markers_10.csv", row.names = TRUE)

#remove cluster 9 and 10, which show cell cycle related genes only
subset <- subset(subset, idents = c("0","1", "2", "3", "4","5", "7", "8"))

#renumber clusters
current.cluster.ids <- c(0,1,2,3,4,5,7,8)

new.cluster.ids <- c("0", "1", "2", "7", "3", "5", "6", "4")

subset@active.ident <- plyr::mapvalues(x = subset@active.ident, from = current.cluster.ids, to = new.cluster.ids)

DimPlot(subset, reduction = "umap", label = TRUE, repel = TRUE)

my_levels <- c("0", "1", "2", "3", "4", "5", "6", "7")

subset@active.ident <- factor(x = subset@active.ident, levels = my_levels)

#figure S7A and S7B
cols<-c("#9E0142", "#084081", "#2ca25f", "#fd8d3c" , "#8c96c6", "#41b6c4", "#525252", "#810f7c" )
DimPlot(subset, reduction = "umap", label = TRUE, cols = cols, pt.size = 1.5)
DimPlot(subset, group.by = "orig.ident", label = TRUE, repel = TRUE)

#figure S7C
FeaturePlot(subset, features = c("TNNT2", "ACTN2"), cols = rev(brewer.pal(n = 11, name = "RdBu")))

#figure 7B
all.markers <- FindAllMarkers(object = subset, only.pos = TRUE, min.pct = 0.25,  thresh.use = 0.25)

top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

DoHeatmap(object = subset, features  = top10$gene, label = T, group.colors = cols) + scale_fill_gradientn(colors =rev(brewer.pal(n = 11, name = "RdBu")))


#remove non-CM clusters according to TNNT2 and ACTN2 expression
subset <- subset(subset, idents = c("0","1", "2","3", "4"))

#Fig 7A 
cols<-c("#9E0142", "#084081", "#2ca25f", "#fd8d3c" , "#8c96c6" )
DimPlot(subset, group.by = "orig.ident", label = TRUE, repel = TRUE) + NoLegend()
DimPlot(subset, reduction = "umap", label = TRUE, repel = TRUE, cols=cols)

#Fig. 7D
VlnPlot(object = subset, features =c("NKX2-5", "NPPA", "HAMP", "CPNE5"),cols = cols)

#Fig. S7E
FeaturePlot(object = subset, features =c("PDGFD", "TXLNB", "GPRIN3"))


##highligh SAN-subpopulation identified in Fig.3 as shown in Fig.7D
#extract SANCM differentiation
Idents(subset)<-"orig.ident"

SAN<-subset(subset, idents= c("D19SAN_1", "D19SAN_2", "D19SAN_3"))
DimPlot(SAN, reduction = "umap", label = TRUE, pt.size = 0.75)

#Normalize data
data.list <- SplitObject(SAN, split.by = "batch")

##Prior to finding anchors, we perform standard preprocessing, 
##and identify variable features individually for each. Note that Seurat v3 
##implements an improved method for variable feature selection based on a 
##variance stabilizing transformation ("vst")

for (i in 1:length(data.list)) { data.list[[i]] <- SCTransform(data.list[[i]], verbose = TRUE)}

#select features for downstream integration, and run PrepSCTIntegration, which 
#ensures that all necessary Pearson residuals have been calculated

data.features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = 3000)
data.list <- PrepSCTIntegration(object.list = data.list, anchor.features = data.features, 
                                verbose = TRUE)

#identify anchors and integrate the datasets
#make sure to set normalization.method = 'SCT'

data.anchors <- FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT", 
                                       anchor.features = data.features, verbose = TRUE)
data.integrated <- IntegrateData(anchorset = data.anchors, normalization.method = "SCT",verbose = TRUE)



##we continue working in the assay "integrated"
DefaultAssay(data.integrated) <- "integrated"

#proceed with downstream analysis (i.e. visualization, clustering) on the integrated dataset
#do not scale integrated data

data.integrated <- RunPCA(data.integrated, verbose = FALSE)
DimPlot(data.integrated, reduction = "pca", dims = c(1,2))
DimHeatmap(data.integrated, dims = 1, cells = 500, balanced = TRUE)

ElbowPlot(data.integrated, ndims = 50)

###15 dims were chosed due to p-value in Jackstraw plot
data.integrated <- RunUMAP(data.integrated, dims = 1:20)
DimPlot(data.integrated, reduction = "umap")
data.integrated <- FindNeighbors(data.integrated, dims = 1:20)
data.integrated <- FindClusters(object = data.integrated, reduction = "pca",dims = 1:15, resolution = 0.6,random.seed = 2020)
head(Idents(data.integrated), 5)
data.integrated <- RunUMAP(data.integrated, dims = 1:20)
DimPlot(data.integrated, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(data.integrated, group.by = "orig.ident", label = TRUE, repel = TRUE)

DefaultAssay(data.integrated) <- "RNA"


# Normalize RNA data for visualization purposes
data.integrated <- NormalizeData(data.integrated, verbose = FALSE)
data.integrated <- ScaleData(data.integrated, verbose = FALSE)

VlnPlot(data.integrated, c("SHOX2", "TBX18", "NKX2-5", "NPPA", "KCNIP2", "HAMP", "ADM"))


SANhead <- WhichCells(data.integrated, idents = "0")
SANtail <- WhichCells(data.integrated, idents = "1")
SANTZ <- WhichCells(data.integrated, idents = "2")

DimPlot(subset, label=T, cells.highlight= SANhead, pt.size = 2)
DimPlot(subset, label=T, cells.highlight= SANtail)
DimPlot(subset, label=T, cells.highlight= SANTZ)

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
