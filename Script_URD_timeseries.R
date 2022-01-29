##URD analysis of timerseries dataset (D0-D19 hiPSC-SANCMs)

library(Seurat)
library(URD)
library(Matrix)
library(fossil) 
library(dplyr)
library(plyr)
library(ggplot2)
library(cowplot)
library(patchwork)

#laod seurat object for URD analysis
setwd("/path/to/directory/")
subset<-readRDS("03_subset_20dims.rds")

##prepare seurat object for URD analysis
#cluster numbering was customized
current.cluster.ids <- c(0,1,2,3,4,5,6,8,9,10,11,12)
new.cluster.ids <- c("6", "5", "8", "4", "0", "2", "1",
                     "7", "9", "5", "3", "10")
subset@active.ident <- plyr::mapvalues(x = subset@active.ident, from = current.cluster.ids, to = new.cluster.ids)
#and ordered
my_levels <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
subset@active.ident <- factor(x = subset@active.ident, levels = my_levels)

#calculate tSNE for URD usign same dimensions
subset <- RunTSNE(subset, dims = 1:20)
DimPlot(subset, reduction = "tsne",cols = getPalette(colourCount))
DimPlot(subset, reduction = "tsne", group.by = "orig.ident", label = TRUE, repel = TRUE)
saveRDS(subset, "06_subset_forURD.rds") ##cluster remain numbered

#remove unwanted clusters for URD analysis, such as hiPSC clusters (cl0, cl1) and 
#endodermal-like cluster (cl3), since we reasoned that diversification will not occur before
#mesoderm stage (cl2)
data<-readRDS("06_subset_forURD.rds")
DimPlot(data, reduction = "umap", label = TRUE, repel = TRUE)
data <- SubsetData(object = data, ident.remove = c("0", "1", "3"))
DimPlot(data, reduction = "umap", label = TRUE, repel = TRUE)

#save dataset for URD analysis
setwd("/path/to/directory")
saveRDS(data, "input_URD.rds")

##convert Seurat object to URD object (for seurat v3: https://github.com/farrellja/URD/issues/29)
data<-readRDS("input_URD.rds")

seuratToURD2 <- function(seurat.object) {
  if (requireNamespace("Seurat", quietly = TRUE)) {
    # Create an empty URD object
    ds <- new("URD")
    
    # Copy over data
    ds@logupx.data <- as(as.matrix(seurat.object@assays$RNA@data), "dgCMatrix")
    if(!any(dim(seurat.object@assays$RNA@counts) == 0)) ds@count.data <- as(as.matrix(seurat.object@assays$RNA@counts[rownames(seurat.object@assays$RNA@data), colnames(seurat.object@assays$RNA@data)]), "dgCMatrix")
    
    # Copy over metadata
    ## TO DO - grab kmeans clustering info
    get.data <- NULL
    if (.hasSlot(seurat.object, "data.info")) { 
      get.data <- as.data.frame(seurat.object@assays$RNA@data.info)
    } else if (.hasSlot(seurat.object, "meta.data")) { 
      get.data <- as.data.frame(seurat.object@meta.data) 
    }
    if(!is.null(get.data)) {
      di <- colnames(get.data)
      m <- grep("res|cluster|Res|Cluster", di, value=T, invert = T) # Put as metadata if it's not the result of a clustering.
      discrete <- apply(get.data, 2, function(x) length(unique(x)) / length(x))
      gi <- di[which(discrete <= 0.015)]
      ds@meta <- get.data[,m,drop=F]
      ds@group.ids <- get.data[,gi,drop=F]
    }
    
    # Copy over var.genes
    if(length(seurat.object@assays$integrated@var.features > 0)) ds@var.genes <- seurat.object@assays$integrated@var.features
    
    # Move over tSNE projection
    if (.hasSlot(seurat.object, "tsne.rot")) {
      if(!any(dim(seurat.object@tsne.rot) == 0)) {
        ds@tsne.y <- as.data.frame(seurat.object@tsne.rot)
        colnames(ds@tsne.y) <- c("tSNE1", "tSNE2")
      }
    } else if (.hasSlot(seurat.object, "reductions")) {
      if(("tsne" %in% names(seurat.object@reductions)) && !any(dim(seurat.object@reductions$tsne) == 0)) {
        ds@tsne.y <- as.data.frame(seurat.object@reductions$tsne@cell.embeddings)
        colnames(ds@tsne.y) <- c("tSNE1", "tSNE2")
      }
    }
    
    # Move over PCA results
    if (.hasSlot(seurat.object, "pca.x")) {
      if(!any(dim(seurat.object@pca.x) == 0)) {
        ds@pca.load <- seurat.object@pca.x
        ds@pca.scores <- seurat.object@pca.rot
        warning("Need to set which PCs are significant in @pca.sig")
      }
      ## TO DO: Convert SVD to sdev
    } else if (.hasSlot(seurat.object, "reductions")) {
      if(("pca" %in% names(seurat.object@reductions)) && !any(dim(Loadings(seurat.object, reduction = "pca")) == 0)) {
        ds@pca.load <- as.data.frame(Loadings(seurat.object, reduction = "pca"))
        ds@pca.scores <- as.data.frame(seurat.object@reductions$pca@cell.embeddings)
        ds@pca.sdev <- seurat.object@reductions$pca@stdev
        ds@pca.sig <- pcaMarchenkoPastur(M=dim(ds@pca.scores)[1], N=dim(ds@pca.load)[1], pca.sdev=ds@pca.sdev)
      }
    }
    return(ds)
  } else {
    stop("Package Seurat is required for this function. To install: install.packages('Seurat')\n")
  }
} ## https://github.com/farrellja/URD/issues/29

data<-seuratToURD2(data)

#annotate stages
data@group.ids$STAGE<-data@group.ids$orig.ident
data@meta$STAGE<-data@group.ids$orig.ident

#clean stage annotation
stages<- data@group.ids$STAGE
new.stages<- gsub("_.*","", stages)
table(new.stages)
new.stages<- gsub("D4","D04", new.stages)
new.stages<- gsub("D5","D05", new.stages)
new.stages<- gsub("D5PM","D05", new.stages)
new.stages<- gsub("D6","D06", new.stages)
data@group.ids$STAGE<- new.stages
data@meta$STAGE<- new.stages

##URD analysis
#remove outliers
#identify cells that are poorly connected
# calculate a k-nearest neighbor graph
data <- calcKNN(data, nn=100)

#plot cells according to their distance to their nearest and 20th nearest neighbors, and identify those with unusually large distances.
outliers <- knnOutliers(data, nn.1=1, nn.2=20, x.max=30, slope.r=1.1, int.r=2.9, slope.b=0.85, int.b=10, title = "Identifying Outliers by k-NN Distance.")

#we cropped cells to the right of the green line (those that are unusually far from 
#their nearest neighbor) and cells above the blue or red lines (those that are unusually far 
#from their 20th nearest neighbor, given their distance to their 1st nearest neighbor).
#subset object to eliminate outliers
cells.keep <- setdiff(colnames(data@logupx.data), outliers)
data <- urdSubset(data, cells.keep=cells.keep)

#calculate diffusion map
data
#URD object: 18834 genes x 2432 cells. knn=50
data <- calcDM(data, knn=50, sigma.use= 12)
plotDimArray(object = data, reduction.use = "dm", dims.to.plot = 1:18, label="STAGE", plot.title="", outer.title="STAGE - Diffusion Map", legend=F, alpha=0.3)
plotDim(data, "STAGE")

##Pseudotime
#calculate pseudotime "floods"
#define the root cells as cells at mesodermal stage (D04)
root.cells <- rownames(data@meta)[data@meta$STAGE=="D04"]
#do the flood
flood.result <- floodPseudotime(data, root.cells=root.cells, n=100, minimum.cells.flooded=2, verbose=T)
#then we process the simulations into a pseudotime
data <- floodPseudotimeProcess(data, flood.result, floods.name="pseudotime")
#check whether enough simulations were run
pseudotimePlotStabilityOverall(data)
#inspect pseudotime
plotDim(data, "pseudotime", plot.title = "Pseudotime")
#create a data.frame that includes pseudotime and stage information
gg.data <- cbind(data@pseudotime, data@meta[rownames(data@pseudotime),])
#plot
ggplot(gg.data, aes(x=pseudotime, color=STAGE, fill=STAGE)) + geom_density(alpha=0.4) + theme_bw()


#Clustering and determining tips
library("knitr")
opts_chunk$set(tidy.opts=list(width.cutoff=80),tidy=TRUE, dev="png", dpi=150)
library(URD)
library(gridExtra) # grid.arrange
#Trim to cells from final stage
#perform secondary clustering of D19 SANCM dataset
#import seurat object and subset D19 dataset
setwd("/path/to/directory/")
D19<-readRDS("subset_forURD.rds")
DimPlot(D19, reduction = "tsne",cols = getPalette(colourCount))
Idents(D19)<-"orig.ident"
D19<-subset(D19, idents = c("D19_1", "D19_2","D19_3"))

#normalization using SCTtransform + batch correction
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
#we set normalization.method = 'SCT'
data.anchors <- FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT", 
                                       anchor.features = data.features, verbose = TRUE)
data.integrated <- IntegrateData(anchorset = data.anchors, normalization.method = "SCT",verbose = TRUE)
#we continue working in the assay "integrated"
DefaultAssay(data.integrated) <- "integrated"
#proceed with downstream analysis (i.e. visualization, clustering) on the integrated datasetdata.integrated <- RunPCA(data.integrated, verbose = FALSE)
data.integrated <- RunPCA(data.integrated, verbose = FALSE)
#elbowplot was additionally used to determine dimensionality of dataset
ElbowPlot(data.integrated, ndims = 50)
#the top 20 PC were chosen for clustering
data.integrated <- FindNeighbors(data.integrated, dims = 1:20)
data.integrated <- FindClusters(object = data.integrated, reduction = "pca",dims = 1:20, resolution = 0.6,random.seed = 2020)
head(Idents(data.integrated), 5)
data.integrated <- RunUMAP(data.integrated, dims = 1:20)
DimPlot(data.integrated, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(data.integrated, group.by = "orig.ident", label = TRUE, repel = TRUE)

#for differential gene expression after integration, we switch back to "RNA" assay and perform logNorm
DefaultAssay(data.integrated) <- "RNA"

#normalize RNA data for visualization purposes
data.integrated <- NormalizeData(data.integrated, verbose = FALSE)
data.integrated <- ScaleData(data.integrated, verbose = FALSE)

#clusters enriched in ERCC-spike-in DNA can be identified as cells with very low counts
VlnPlot(data.integrated, "nFeature_RNA")

##cluster 3 showed very low count of genes and was therefore removed for further analysis
subset <- SubsetData(object = data.integrated, ident.remove = "3")

##sava D19 dataset for tip cluster annotation
saveRDS(subset, "subset_D19_SAN_20dims_annot.rds")

##import D19 dataset for tips annotation
D19<-readRDS("subset_D19_SAN_20dims_annot.rds")
#convert D19 seurat object to URD object for tips annotation
D19<-seuratToURD2(D19)

#transfer clustering to main object (in new column D19)
data@group.ids[rownames(D19@group.ids), "D19"] <- D19@group.ids$seurat_clusters
#remove cells which could not be annotated (NA)
tail(data@group.ids)
groupids<-data@group.ids
groupids<-groupids %>% drop_na(seurat_clusters)
tail(groupids)
data@group.ids <- groupids
tail(data@group.ids, 50)
#check whether tips are indeed tips (in this case SANhead cluster)
data@group.ids$pop <- NA
data@group.ids[cellsInCluster(data, "D19", "2"), "pop"] <- "1"
plotDim(data, label = "pop", plot.title = "SANhead", reduction.use = "dm", dim.x = 5, dim.y = 6, legend = TRUE, alpha= 0.35)

#biased random walks
#copy cluster identities from D19 object to a new clustering ("tip.clusters") in the full URD object.
data@group.ids[rownames(D19@group.ids), "tip.clusters"] <- D19@group.ids$seurat_clusters

#determine the parameters of the logistic used to bias the transition probabilities. The procedure
#is relatively robust to this parameter, but the cell numbers may need to be modified for larger
#or smaller data sets.
data.ptlogistic <- pseudotimeDetermineLogistic(data, "pseudotime", optimal.cells.forward=50, max.cells.back=80, do.plot = T)
#bias the transition matrix according to pseudotime
data.biased.tm <- as.matrix(pseudotimeWeightTransitionMatrix(data, "pseudotime", logistic.params=data.ptlogistic))
#simulate the biased random walks from each tip
data.walks <- simulateRandomWalksFromTips(data, tip.group.id="tip.clusters", root.cells=root.cells, transition.matrix = data.biased.tm, n.per.tip = 25000, root.visits = 1, max.steps = 5000, verbose = F)
#process the biased random walks into visitation frequencies
data <- processRandomWalksFromTips(data, data.walks, verbose = F)
#visualize the tips and the visitation of cells from each tip on the dataset
plotDim(data, "tip.clusters", plot.title="Cells in each tip")
plotDim(data, "visitfreq.log.pseudotime", plot.title="Visitation frequency", transitions.plot=10000)

#build tree
#load the cells used for each tip into the URD object
data.tree <- loadTipCells(data, "tip.clusters")

#build the tree using cluster 0, 1, 2, 5 as tips
data.tree <- buildTree(data.tree, pseudotime = "pseudotime", tips.use= c("0", "1", "2", "5"), divergence.method = "preference", cells.per.pseudotime.bin = 35, bins.per.pseudotime.window = 10, save.all.breakpoint.info = T, p.thresh=0.000001)

#name the segments based on our previous determination of the identity of tips 1 and 2.
data.tree <- nameSegments(data.tree, segments=c("0", "1", "2", "5"), segment.names = c("SANTZ", "SANhead", "SANtail", "proepicardial calls"), short.names = c("SANTZ", "SANhead", "SANtail", "proepi"))

#plot meta data or gene expression in dendrogram recovered by URD
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(8, "Spectral"))
colourCount = 5

#pseudotime tree as shown in Fig5A
plotTree(data.tree, "STAGE", title="Developmental Stage",legend.point.size = 5, cell.size = 0.9, cell.alpha = 0.5, tree.size = 1, tree.alpha = 0.7, discrete.colors = getPalette(colourCount))

#visualisation different genes in dendrogramm Fig5B and C
plotTree(data.tree1, "TNNT2", legend.point.size = 8* cell.size, cell.size = 0.9, cell.alpha = 0.5, tree.size = 1.5)
plotTree(data.tree1, "MYH6", legend.point.size = 8* cell.size, cell.size = 0.9, cell.alpha = 0.5, tree.size = 1.5)
plotTree(data.tree1, "WT1", legend.point.size = 8* cell.size, cell.size = 0.9, cell.alpha = 0.5, tree.size = 1.5)
plotTree(data.tree1, "ALDH1A2", legend.point.size = 8* cell.size, cell.size = 0.9, cell.alpha = 0.5, tree.size = 1.5)

#visualisation different genes in dendrogramm Fig6A and C
plotTree(data.tree1, "BMP4", legend.point.size = 8* cell.size, cell.size = 0.9, cell.alpha = 0.5, tree.size = 1.5)
plotTree(data.tree1, "BMP2", legend.point.size = 8* cell.size, cell.size = 0.9, cell.alpha = 0.5, tree.size = 1.5)
plotTree(data.tree1, "HTRA1", legend.point.size = 8* cell.size, cell.size = 0.9, cell.alpha = 0.5, tree.size = 1.5)
plotTree(data.tree1, "FBN2", legend.point.size = 8* cell.size, cell.size = 0.9, cell.alpha = 0.5, tree.size = 1.5)
plotTree(data.tree1, "DKK1", legend.point.size = 8* cell.size, cell.size = 0.9, cell.alpha = 0.5, tree.size = 1.5)
plotTree(data.tree1, "WNT5A", legend.point.size = 8* cell.size, cell.size = 0.9, cell.alpha = 0.5, tree.size = 1.5)
plotTree(data.tree1, "SFRP1", legend.point.size = 8* cell.size, cell.size = 0.9, cell.alpha = 0.5, tree.size = 1.5)
plotTree(data.tree1, "APP", legend.point.size = 8* cell.size, cell.size = 0.9, cell.alpha = 0.5, tree.size = 1.5)


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

