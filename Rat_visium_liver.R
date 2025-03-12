library(Seurat)
devtools::install_github('satijalab/seurat-data')
library(SeuratData)
library(dplyr)
library(ggplot2)
library(patchwork)

visium_ct <- Load10X_Spatial("/Users/marina/Documents/Projectes/Visium 10x/Rat_visium_liver/spaceranger/Control/", 
                          filename = "filtered_feature_bc_matrix.h5", 
                          assay = "Spatial", 
                          slice = "slice1",
                          filter.matrix = TRUE)
visium_ch <- Load10X_Spatial("/Users/marina/Documents/Projectes/Visium 10x/Rat_visium_liver/spaceranger/Cirrhotic/", 
                          filename = "filtered_feature_bc_matrix.h5", 
                          assay = "Spatial", 
                          slice = "slice2",
                          filter.matrix = TRUE)
visium_r_1<- Load10X_Spatial("/Users/marina/Documents/Projectes/Visium 10x/Rat_visium_liver/spaceranger/Regression/R_1/", 
                        filename = "filtered_feature_bc_matrix.h5", 
                        assay = "Spatial", 
                        slice = "slice3",
                        filter.matrix = TRUE)
visium_r_2 <- Load10X_Spatial("/Users/marina/Documents/Projectes/Visium 10x/Rat_visium_liver/spaceranger/Regression/R_2/", 
                        filename = "filtered_feature_bc_matrix.h5", 
                        assay = "Spatial", 
                        slice = "slice4",
                        filter.matrix = TRUE)

#Samples name
visium_ct$sample <-"Control"
visium_ch$sample <-"Cirrhotic"
visium_r_1$sample <-"Regression_1"
visium_r_2$sample <-"Regression_2"

#MITO

merged_visium <- SetIdent(merged_visium, value = merged_visium@meta.data$seurat_clusters)

visium_ct$mitoRatio <- PercentageFeatureSet(object = visium_ct, pattern = "Mt-")
visium_ct$mitoRatio <- visium_ct@meta.data$mitoRatio/100
SpatialFeaturePlot(visium_ct, features = "mitoRatio", alpha = c(0.1,1))
VlnPlot(visium_ct, features = "mitoRatio")

visium_ch$mitoRatio <- PercentageFeatureSet(object = visium_ch, pattern = "Mt-")
visium_ch$mitoRatio <- visium_ch@meta.data$mitoRatio/100
SpatialFeaturePlot(visium_ch, features = "mitoRatio", alpha = c(0.1,1))
VlnPlot(visium_ch, features = "mitoRatio")

visium_r_1$mitoRatio <- PercentageFeatureSet(object = visium_r_1, pattern = "Mt-")
visium_r_1$mitoRatio <- visium_r_1@meta.data$mitoRatio/100
SpatialFeaturePlot(visium_r_1, features = "mitoRatio", alpha = c(0.1,1))
VlnPlot(visium_r_1, features = "mitoRatio")

visium_r_2$mitoRatio <- PercentageFeatureSet(object = visium_r_2, pattern = "Mt-")
visium_r_2$mitoRatio <- visium_r_2@meta.data$mitoRatio/100
SpatialFeaturePlot(visium_r_2, features = "mitoRatio", alpha = c(0.1,1))
VlnPlot(visium_r_2, features = "mitoRatio")

#Split Seurat object 
DefaultAssay(all.liver.combined) <- "RNA"
all.liver.combined <- SetIdent(all.liver.combined, value = all.liver.combined@meta.data$orig.ident)
control <- subset(all.liver.combined, idents = c("Control"))
cirrhotic <- subset(all.liver.combined, idents = c("Cirrhotic"))
regression <- subset(all.liver.combined, idents = c("Regression"))

#MERGE
merged_visium <-merge(visium_ct, c(visium_ch,visium_r_1,visium_r_2), add.cell.id = c("Control", "Cirrhotic", "Regression_1", "Regression_2"))
SpatialFeaturePlot(merged_visium, features = "nFeature_Spatial")
VlnPlot(merged_visium, features = c("nCount_Spatial", "nFeature_Spatial", "mitoRatio"), split.by = "sample", pt.size = 0)
merged_visium<- SetIdent(merged_visium, value = merged_visium@meta.data$sample)


#SCTransform builds regularized negative binomial models of gene expression in order 
#to account for technical artifacts while preserving biological variance. NORMALIZATION

visium_ct <- SCTransform(visium_ct, assay = "Spatial", verbose = FALSE)
visium_ch <- SCTransform(visium_ch, assay = "Spatial", verbose = FALSE)
visium_r_1 <- SCTransform(visium_r_1, assay = "Spatial", verbose = FALSE)
visium_r_2 <- SCTransform(visium_r_2, assay = "Spatial", verbose = FALSE)
merged_visium <- SCTransform(merged_visium, assay = "Spatial", verbose = FALSE)

#Gene expression Visualization

SpatialFeaturePlot(merged_visium, c("Apoe", "Ambp"), alpha = c(0.1,1))
SpatialFeaturePlot(visium_ch, c("Glul", "Hal", "Alb", "Cd74"))
SpatialFeaturePlot(visium_r_1, c("C3", "Hp", "Scd", "App"))
SpatialFeaturePlot(visium_r_2, c("C3", "Hp", "Scd", "App"))

#Dimensionality reduction, clustering and visualization
merged_visium <- RunPCA(merged_visium, assay = "SCT", verbose = FALSE)
merged_visium <- FindNeighbors(merged_visium, reduction = "pca", dims = 1:30)
merged_visium <- FindClusters(merged_visium, verbose = FALSE, resolution = 1.5)
merged_visium <- RunUMAP(merged_visium, reduction = "pca", dims = 1:30)

DimPlot(merged_visium, reduction = "umap", label = T, group.by = "seurat_clusters")
SpatialDimPlot(merged_visium, label.size = 3, alpha = c(0.5,1))
SpatialDimPlot(visium_r_1)
VlnPlot(visium_r_1, features = c("Apoa1", "Pnpla3", "Scd", "App"), ncol = 2)

visium_ct <- RunPCA(visium_ct, assay = "SCT", verbose = FALSE)
visium_ct <- FindNeighbors(visium_ct, reduction = "pca", dims = 1:20)
visium_ct <- FindClusters(visium_ct, verbose = FALSE, resolution = 1.5)
visium_ct <- RunUMAP(visium_ct, reduction = "pca", dims = 1:30)

DimPlot(visium_ct, reduction = "umap", label = T)
SpatialDimPlot(visium_ct, label.size = 3, alpha = c(0.3,1))
VlnPlot(visium_ct, features = c("Cd74"))
ct_markers <- FindAllMarkers(visium_ct, min.pct = 0.25)

visium_ch <- RunPCA(visium_ch, assay = "SCT", verbose = FALSE)
visium_ch <- FindNeighbors(visium_ch, reduction = "pca", dims = 1:10)
visium_ch <- FindClusters(visium_ch, verbose = FALSE, resolution = 1.5)
visium_ch <- RunUMAP(visium_ch, reduction = "pca", dims = 1:30)

DimPlot(visium_ch, reduction = "umap", label = T)
SpatialDimPlot(visium_ch, label.size = 3, alpha = c(0.3,1), label = T)
VlnPlot(visium_ch, features = "nCount_Spatial")
SpatialDimPlot(visium_ch, cells.highlight = CellsByIdentities(object = visium_ch, idents = c(2, 1, 4, 3, 5, 8)), facet.highlight = TRUE, ncol = 3)
ch_markers <- FindAllMarkers(visium_ch, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 2)

visium_r_1 <- RunPCA(visium_r_1, assay = "SCT", verbose = FALSE)
visium_r_1 <- FindNeighbors(visium_r_1, reduction = "pca", dims = 1:30)
visium_r_1 <- FindClusters(visium_r_1, verbose = FALSE, resolution = 1.5)
visium_r_1 <- RunUMAP(visium_r_1, reduction = "pca", dims = 1:30)

DimPlot(visium_r_1, reduction = "umap", label = T)
SpatialDimPlot(visium_r_1, label.size = 3, alpha = c(0.3,1))
SpatialDimPlot(visium_r_1)
VlnPlot(visium_r_1, features = c("Apoa1", "Pnpla3", "Scd", "App"), ncol = 2)

visium_r_2 <- RunPCA(visium_r_2, assay = "SCT", verbose = FALSE)
visium_r_2 <- FindNeighbors(visium_r_2, reduction = "pca", dims = 1:30)
visium_r_2 <- FindClusters(visium_r_2, verbose = FALSE, resolution = 1.5)
visium_r_2 <- RunUMAP(visium_r_2, reduction = "pca", dims = 1:30)

DimPlot(visium_r_2, reduction = "umap", label = T)
SpatialDimPlot(visium_r_2, label.size = 3, alpha = c(0.3,1))
SpatialDimPlot(visium_r_2)
VlnPlot(visium_r_2, features = "nCount_Spatial")

#INTEGRATON from https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_07_spatial.html
#create a list of the original data that we loaded
st.list = list(Control=visium_ct, Cirrhotic= visium_ch, Regression_1 = visium_r_1, Regression_2= visium_r_2)
#run SCT on all datasets

st.list = lapply(st.list, SCTransform, assay= "Spatial", method = "poisson")
options(future.globals.maxSize = 2000 * 1024^2)
st.features = SelectIntegrationFeatures(st.list, nfeatures = 3000, verbose = FALSE)
st.list <- PrepSCTIntegration(object.list = st.list, anchor.features = st.features, verbose = FALSE)
#perform integration
int.anchors <- FindIntegrationAnchors(object.list = st.list, normalization.method = "SCT", verbose = FALSE, anchor.features = st.features)
visium.integrated <- IntegrateData(anchorset = int.anchors, normalization.method = "SCT", verbose = FALSE)
rm(int.anchors, st.list)

#Run dimensionality reduction and clustering as before
visium.integrated <- RunPCA(visium.integrated, verbose = FALSE)
visium.integrated <- FindNeighbors(visium.integrated, dims = 1:30)
visium.integrated <- FindClusters(visium.integrated, verbose = FALSE)
visium.integrated <- RunUMAP(visium.integrated, dims = 1:30)
DimPlot(visium.integrated, reduction = "umap", group.by = c("ident", "sample"))
SpatialDimPlot(visium.integrated)

DimPlot(visium.integrated, reduction = "umap", label = T, split.by = "sample")
SpatialDimPlot(visium.integrated, label.size = 3, ncol = 2, alpha = c(0.3,1), facet.highlight = TRUE) 
SpatialDimPlot(visium.integrated, cells.highlight = CellsByIdentities(object = visium.integrated, idents = c(9)), ncol = 2)

visium.integrated$sample <-factor(x = visium.integrated$sample, levels = c("Control", "Cirrhotic", "Regression_1", "Regression_2")) #Reordenar

VlnPlot(visium.integrated, features = c("Pdgfrb", "Cd68", "Cd163", "Mrc1"), ncol = 2, pt.size = 0)
VlnPlot(visium.integrated, features = c("nCount_Spatial", "nFeature_Spatial", "mitoRatio"), split.by = "seurat_clusters", pt.size = 0)

#SUbset without cluster 2 and 9
visium.integrated_subset <- subset(x = visium.integrated, idents = c("0","1","3","4","5","6","7","8","10","11","12"))


#RENAME
new_clusters_vis <- c("HSCs", "Periportal Hepatocytes", "", "Erytroid Cells", "Pericentral Hepatocytes", "Interzonal Hepatocytes", "Monocytes/LSEC?", "B-cell", "Erytroid cells", "Hepatocytes",  )
visium.integrated <- RenameIdents(visium.integrated, new_clusters_vis)
visium.integrated <- SetIdent(visium.integrated, value = visium.integrated@meta.data$new_clusters_vis)

#Identification of spatially variable features

# differential expression 
DefaultAssay(visium.integrated) <- "integrated"
int_markers_subset <- FindAllMarkers(visium.integrated_subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 2)
int_de_12 <- FindMarkers(visium.integrated, ident.1 = 12, ident.2 = NULL, only.pos = T, min.pct = 0.25, logfc.threshold = 2, p)

#ordenar per FoldChange
int_de_0 <-int_de_0[order(-int_de_0$avg_log2FC),]

# plot top markers
SpatialFeaturePlot(object = visium_ch, features = rownames(de_ch_8)[1:6],
                   alpha = c(0.1, 1), ncol = 3)

de_visium_ct <- FindSpatiallyVariableFeatures(visium_ct, assay = "SCT", features = VariableFeatures(visium_ct)[1:1000],selection.method = "markvariogram")
top.features_ct <- head(SpatiallyVariableFeatures(de_visium_ct, selection.method = "markvariogram"),6)
SpatialFeaturePlot(de_visium_ct, features = c("Cd74", "Ccl5"), alpha = c(0.1,1))

de_visium_ch <- FindSpatiallyVariableFeatures(visium_ch, assay = "SCT", features = VariableFeatures(visium_ch)[1:1000],selection.method = "markvariogram")
top.features_ch <- head(SpatiallyVariableFeatures(de_visium_ch, selection.method = "markvariogram"),6)
SpatialFeaturePlot(de_visium_ch, features = c("Cd74", "Ccl5"), alpha = c(0.1,1))

de_visium_r1 <- FindSpatiallyVariableFeatures(visium_r_1, assay = "SCT", features = VariableFeatures(visium_r_1)[1:1000],selection.method = "markvariogram")
top.features_r1 <- head(SpatiallyVariableFeatures(de_visium_r1, selection.method = "markvariogram"),6)
SpatialFeaturePlot(de_visium_r1, features = c("Cd74", "Ccl5"), alpha = c(0.1,1))

de_visium_r2 <- FindSpatiallyVariableFeatures(visium_r_2, assay = "SCT", features = VariableFeatures(visium_r_2)[1:1000],selection.method = "markvariogram")
top.features_r2 <- head(SpatiallyVariableFeatures(de_visium_r2, selection.method = "markvariogram"),10)
SpatialFeaturePlot(de_visium_r2, features = top.features_r2, alpha = c(0.1,1))
SpatialFeaturePlot(de_visium_r2, features = c("Cd74", "Ccl5"), alpha = c(0.1,1))

dif_visium_ct<- FindAllMarkers(visium_ct, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
dif_visium_ct.list <- dif_visium_ct%>%group_by(cluster) %>% top_n(n=10, wt = avg_log2FC) %>% select(cluster, marker_gene = gene)

dif_visium_ch<- FindAllMarkers(visium_ch, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
dif_visium_ch.list <- dif_visium_ch%>%group_by(cluster) %>% top_n(n=10, wt = avg_log2FC) %>% select(cluster, marker_gene = gene)

dif_visium_r1<- FindAllMarkers(visium_r_1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
dif_visium_r1.list <- dif_visium_r1%>%group_by(cluster) %>% top_n(n=10, wt = avg_log2FC) %>% select(cluster, marker_gene = gene)

dif_visium_r2<- FindAllMarkers(visium_r_2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
dif_visium_r2.list <- dif_visium_r2%>%group_by(cluster) %>% top_n(n=10, wt = avg_log2FC) %>% select(cluster, marker_gene = gene)

#LABEL TRANSFER FROM SEURAT

DefaultAssay(visium.integrated) <- "integrated"
DimPlot(all.liver.combined, group.by = "new_clusters_id")
all.liver.combined <- SCTransform(all.liver.combined, ncells = 3000, verbose = FALSE) %>% RunPCA(verbose = FALSE) %>% RunUMAP(dims=1:30)

###Potser probar all.liver.combined per a label transfer per a cada dataset individual
DimPlot(control, group.by = "new_clusters_id")
DimPlot(cirrhotic, group.by = "new_clusters_id")
DimPlot(regression, group.by = "new_clusters_id")
control <- SCTransform(control, ncells = 3000, verbose = FALSE) %>% RunPCA(verbose = FALSE) %>% RunUMAP(dims=1:30)

#integrated dataset
anchors_int <- FindTransferAnchors(reference = all.liver.combined, query = visium.integrated, normalization.method = "SCT")
predictions.assay_int <- TransferData(anchorset = anchors_int, refdata = all.liver.combined$new_clusters_id, prediction.assay = TRUE, weight.reduction = visium.integrated[["pca"]], dims = 1:30)
visium.integrated[["predictions"]] <- predictions.assay_int
DefaultAssay(visium.integrated) <- "predictions"
DefaultAssay(merged_visium) <- "Spatial"
SpatialFeaturePlot(visium.integrated, features = "HSCs", alpha = c(0,1))
SpatialFeaturePlot(visium.integrated, features = "Pericentral Hepatocytes")
SpatialFeaturePlot(visium.integrated, features = "Periportal Hepatocytes")
SpatialFeaturePlot(visium.integrated, features = "LSEC", alpha = c(0.3,1))
SpatialFeaturePlot(visium.integrated, features = "Midportal Hepatocytes")
SpatialFeaturePlot(visium.integrated, features = "Kupffer cells M1")
SpatialFeaturePlot(visium.integrated, features = "Kupffer cells M2")
SpatialFeaturePlot(visium.integrated, features = "Proliferating Hepatocytes")
SpatialFeaturePlot(visium.integrated, features = "T/NK-cells")
SpatialFeaturePlot(visium.integrated, features = "Cholangiocytes")
SpatialFeaturePlot(visium.integrated, features = "SAMs")
SpatialFeaturePlot(visium.integrated, features = c("nCount_Spatial"))
VlnPlot(visium.integrated, features = c("nCount_Spatial"))
DimPlot(visium.integrated)
FeaturePlot(visium.integrated, features = c("SAMs"))



#Filter out per mitoRatio >0.1

#individual data
all.liver.combined <-SetIdent(all.liver.combined, value = all.liver.combined@meta.data$new_clusters_id)
anchors_ct <- FindTransferAnchors(reference = control, query = visium_ct, normalization.method = "SCT")
predictions.assay_ct<- TransferData(anchorset = anchors_ct, refdata = control$new_clusters_id, prediction.assay = TRUE, weight.reduction = visium_ct[["pca"]], dims = 1:30)
visium_ct[["predictions"]] <- predictions.assay_ct
DefaultAssay(visium_ct) <- "predictions"
SpatialFeaturePlot(visium_ct, features = "LSEC", pt.size.factor = 1, alpha = c(0.1,1)) #LSEC, Periportal Hepatocytes, Midportal Hepatocytes, Kupffer cells M1, Kupffer cells M2, 
SpatialFeaturePlot(visium_ct, features = "HSCs", pt.size.factor = 1, alpha = c(0.1,1))
visium_ct <- FindSpatiallyVariableFeatures(visium_ct, assay = "predictions", selection.method = "markvariogram", features = rownames(visium_ct), r.metric = 5, slot = "data")
top.clusters_ct <- head(SpatiallyVariableFeatures(visium_ct), 4)
SpatialPlot(object = visium_ct, features = top.clusters_ct, ncol = 2)

anchors_ch <- FindTransferAnchors(reference = cirrhotic, query = visium_ch, normalization.method = "SCT")
predictions.assay_ch<- TransferData(anchorset = anchors_ch, refdata = all.liver.combined$new_clusters_id, prediction.assay = TRUE, weight.reduction = visium_ch[["pca"]], dims = 1:30)
visium_ch[["predictions"]] <- predictions.assay_ch
DefaultAssay(visium_ch) <- "predictions"
p2 <-SpatialFeaturePlot(visium_ch, features = "LSEC", pt.size.factor = 1, alpha = c(0.1,1))
visium_ch <- FindSpatiallyVariableFeatures(visium_ch, assay = "predictions", selection.method = "markvariogram", features = rownames(visium_ch), r.metric = 5, slot = "data")
top.clusters_ch <- head(SpatiallyVariableFeatures(visium_ch), 4)
SpatialPlot(object = visium_ch, features = top.clusters_ch, ncol = 2)

anchors_r1 <- FindTransferAnchors(reference = all.liver.combined, query = visium_r_1, normalization.method = "SCT")
predictions.assay_r1<- TransferData(anchorset = anchors_r1, refdata = all.liver.combined$new_clusters_id, prediction.assay = TRUE, weight.reduction = visium_r_1[["pca"]], dims = 1:30)
visium_r_1[["predictions"]] <- predictions.assay_r1
DefaultAssay(visium_r_1) <- "predictions"
p3 <-SpatialFeaturePlot(visium_r_1, features = "LSEC", pt.size.factor = 1, alpha = c(0.1,1))
visium_r_1 <- FindSpatiallyVariableFeatures(visium_r_1, assay = "predictions", selection.method = "markvariogram", features = rownames(visium_r_1), r.metric = 5, slot = "data")
top.clusters_r1 <- head(SpatiallyVariableFeatures(visium_r_1), 4)
SpatialPlot(object = visium_r_1, features = top.clusters_r1, ncol = 2)

anchors_r2 <- FindTransferAnchors(reference = all.liver.combined, query = visium_r_2, normalization.method = "SCT")
predictions.assay_r2<- TransferData(anchorset = anchors_r2, refdata = all.liver.combined$new_clusters_id, prediction.assay = TRUE, weight.reduction = visium_r_2[["pca"]], dims = 1:30)
visium_r_2[["predictions"]] <- predictions.assay_r2
DefaultAssay(visium_r_2) <- "predictions"
p4 <-SpatialFeaturePlot(visium_r_2, features = "LSEC", pt.size.factor = 1, alpha = c(0.1,1))
visium_r_2 <- FindSpatiallyVariableFeatures(visium_r_2, assay = "predictions", selection.method = "markvariogram", features = rownames(visium_r_2), r.metric = 5, slot = "data")
top.clusters_r2 <- head(SpatiallyVariableFeatures(visium_r_2), 4)
SpatialPlot(object = visium_r_2, features = top.clusters_r2, ncol = 2)

p1 + p2 + p3 +p4

#Deconvolution: method to estimate the abundance or proportion of different celltypes in a bulkRNAseq dataset using single cell reference.
#Using Using SCDC for deconvolution.
inst = installed.packages()

if (!("xbioc" %in% rownames(inst))) {
  remotes::install_github("renozao/xbioc", dependencies = FALSE)
}
install.packages('L1pack')
if (!("SCDC" %in% rownames(inst))) {
  remotes::install_github("meichendong/SCDC", dependencies = FALSE)
}
#From single cell all.liver.combined dataset select 200 cells per subclass
#set subclass as active.ident
all.liver.combined <-subset(all.liver.combined, cells=WhichCells(all.liver.combined, downsample = 200))

all.liver.combined@active.assay = "RNA"
markers_sc <-FindAllMarkers(all.liver.combined, only.pos = TRUE, logfc.threshold = 0.1,test.use = "wilcox", min.pct = 0.05,
                            min.diff.pct = 0.1, max.cells.per.ident = 200, return.thresh = 0.05, assay = "RNA")

#Filter for genes that are also in the ST data
markers_sc <-markers_sc[markers_sc$gene %in% rownames(merged_visium),]

#Select top 20 genes per cluster, select top fist p-value, then absolute diff in pct, then quota of pct
markers_sc$pct.diff <-markers_sc$pct.1 - markers_sc$pct.2
markers_sc$log.pct.diff <- log2((markers_sc$pct.1 * 99 +1)/(markers_sc$pct.2 * 99 +1))
markers_sc %>%
  group_by(cluster) %>%
  top_n(-100, p_val) %>%
  top_n(50, pct.diff) %>%
  top_n(20, log.pct.diff) ->top20

m_feats <-unique(as.character(top20$gene))

#Create Expression sets
eset_SC <- ExpressionSet(assayData = as.matrix(all.liver.combined@assays$RNA@counts[m_feats,
                        ]), phenoData = AnnotatedDataFrame(all.liver.combined@meta.data))
eset_ST <- ExpressionSet(assayData = as.matrix(merged_visium@assays$Spatial@counts[m_feats,]),phenoData = AnnotatedDataFrame(merged_visium@meta.data))

#Deconvolve
library("SCDC")
deconvolution_crc <- SCDC::SCDC_prop(bulk.eset = eset_ST, sc.eset = eset_SC, ct.varname = "new_clusters_id", ct.sub = as.character(unique(eset_SC$new_clusters_id)))


allen_reference <- readRDS("../data/allen_cortex.rds")

#Label transfer example
# note that setting ncells=3000 normalizes the full dataset but learns noise models on 3k
# cells this speeds up SCTransform dramatically with no loss in performance
library(dplyr)
allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)
# After subsetting, we renormalize cortex
cortex <- SCTransform(cortex, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(verbose = FALSE)
# the annotation is stored in the 'subclass' column of object metadata
DimPlot(allen_reference, group.by = "subclass", label = TRUE)

anchors <- FindTransferAnchors(reference = allen_reference, query = cortex, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$subclass, prediction.assay = TRUE,
                                  weight.reduction = cortex[["pca"]], dims = 1:30)
cortex[["predictions"]] <- predictions.assay

DefaultAssay(cortex) <- "predictions"
SpatialFeaturePlot(cortex, features = c("L2/3 IT", "L4"), pt.size.factor = 1.6, ncol = 2, crop = TRUE)

