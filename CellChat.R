#CellChat

#tutorial 
#load(url("https://ndownloader.figshare.com/files/25950872"))
#load("/Users/jinsuoqin/Documents/CellChat/tutorial/data_humanSkin_CellChat.rda")

devtools::install_github("sqjin/CellChat")
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
install.packages("extrafont")
library(extrafont)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

#Extract the cellchat input files from a Seurat V3 object

dplyr::glimpse(CellChatDB.rat$interaction)
CellChatDB.use <- CellChatDB.rat

all.liver.combined <- SetIdent(all.liver.combined, value = all.liver.combined@meta.data$orig.ident)
data.input <- as.data.frame(all.liver.combined@assays$RNA@data)
cells <- all.liver.combined@meta.data
cells <-cells[Cells(all.liver.combined),]
data.input <- data.input[,row.names(cells)]
meta = all.liver.combined@meta.data
cellchat <- createCellChat(object = as.matrix(data.input), meta = meta, group.by = "new_clusters_id")
cellchat@DB <- CellChatDB.use
cellchat <-subsetData(cellchat)
future::plan("multiprocess", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#optional cellchat2 <- projectData(cellchat2, PPI.human)
cellchat <- computeCommunProb(cellchat)
cellchat_subset <- filterCommunication(cellchat, min.cells = 10)
cellchat_subset <- computeCommunProbPathway(cellchat_subset)
cellchat_subset <- aggregateNet(cellchat_subset)
groupSize <- as.numeric(table(cellchat_subset@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_subset@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_subset@net$weight, vertex.weight = groupSize_ct, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

#Heatmap
pathways.show <- c("COLLAGEN")
par(mfrow=c(1,1))
netVisual_heatmap(cellchat_subset, signaling = pathways.show, color.heatmap = "Reds")

#contribution
cellchat_subset <- netAnalysis_computeCentrality(cellchat_subset)
netAnalysis_contribution(cellchat_subset, signaling = pathways.show)
netAnalysis_signalingRole_heatmap(cellchat_subset, signaling = NULL, color.heatmap = "Reds")

plotGeneExpression(cellchat_subset, signaling = "COLLAGEN", enriched.only = FALSE)

#SEPARATING BY CONDITION
control <- subset(x = all.liver.combined, idents = c("Control"))
cirrhotic <- subset(x = all.liver.combined, idents = c("Cirrhotic"))
regression <- subset(x = all.liver.combined, idents = c("Regression"))

#control
control <- SetIdent(control, value = control@meta.data$new_clusters_id)
data.input_ct <- as.data.frame(control@assays$RNA@data)
cells_ct <- control@meta.data
cells_ct <-cells_ct[Cells(control),]
data.input_ct <- data.input_ct[,row.names(cells_ct)]
meta_ct = control@meta.data
cellchat_ct <- createCellChat(object = as.matrix(data.input_ct), meta = meta_ct, group.by = "new_clusters_id")
cellchat_ct@DB <- CellChatDB.use
cellchat_ct <-subsetData(cellchat_ct)
cellchat_ct <- identifyOverExpressedGenes(cellchat_ct)
cellchat_ct <- identifyOverExpressedInteractions(cellchat_ct)
cellchat_ct <- computeCommunProb(cellchat_ct)
cellchat_subset_ct <- filterCommunication(cellchat_ct, min.cells = 10)
cellchat_subset_ct <- computeCommunProbPathway(cellchat_subset_ct)
cellchat_subset_ct <- aggregateNet(cellchat_subset_ct)
groupSize_ct <- as.numeric(table(cellchat_subset_ct@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_subset_ct@net$count, vertex.weight = groupSize_ct, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_subset_ct@net$weight, vertex.weight = groupSize_ct, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
#Individuals
mat <- cellchat_subset_ct@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize_ct, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
#Heatmap
#pathways = MHC, MHC-I, CLEC, MIF, CD99, GALECTIN, APP, ITGB2, LCK, IL16, ANNEXIN, CCL, ADGRE5, CD45, CXCL
#XCR, IFN-II, CD22, TNF, SELPLG, COLLAGEN, PECAM1, MK, BAFF, ICAM, PARs, FN1, BTLA, LAMININ, VISFATIN, SELL
#CD23, CD40, COMPLEMENT, LIGHT, CADM, ESAM, JAM, GAS, CD86, FLT3,VCAM, ALCAM,CD6, VEGF, CHEMERIN, CDH5, EPHB, SEMA4, NCAM
pathways.show <- c("COLLAGEN")
par(mfrow=c(1,1))
netVisual_heatmap(cellchat_subset_ct, signaling = pathways.show, color.heatmap = "Reds")

#contribution
cellchat_subset_ct <- netAnalysis_computeCentrality(cellchat_subset_ct)
netAnalysis_signalingRole_heatmap(cellchat_subset_ct, signaling = NULL, color.heatmap = "Reds", font.size = 6)

netAnalysis_contribution(cellchat_subset_ct, signaling = pathways.show)

#plotgeneExpression
plotGeneExpression(cellchat_subset, signaling = "COLLAGEN", enriched.only = FALSE)

# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat_subset_ct, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
#> [[1]]
# Circle plot
netVisual_individual(cellchat_subset, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
###########
#cirrhotic#
###########
cirrhotic <- SetIdent(cirrhotic, value = cirrhotic@meta.data$new_clusters_id)
data.input_ch <- as.data.frame(cirrhotic@assays$RNA@data)
cells_ch <- cirrhotic@meta.data
cells_ch <-cells_ch[Cells(cirrhotic),]
data.input_ch <- data.input_ch[,row.names(cells_ch)]
meta_ch = cirrhotic@meta.data
cellchat_ch <- createCellChat(object = as.matrix(data.input_ch), meta = meta_ch, group.by = "new_clusters_id")
cellchat_ch@DB <- CellChatDB.use
cellchat_ch <-subsetData(cellchat_ch)
cellchat_ch <- identifyOverExpressedGenes(cellchat_ch)
cellchat_ch <- identifyOverExpressedInteractions(cellchat_ch)
cellchat_ch <- computeCommunProb(cellchat_ch)
cellchat_subset_ch <- filterCommunication(cellchat_ch, min.cells = 10)
cellchat_subset_ch <- computeCommunProbPathway(cellchat_subset_ch)
cellchat_subset_ch <- aggregateNet(cellchat_subset_ch)

groupSize_ch <- as.numeric(table(cellchat_subset_ch@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_subset_ch@net$count, vertex.weight = groupSize_ch, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_subset_ch@net$weight, vertex.weight = groupSize_ch, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

#Individuals
mat <- cellchat_subset_ch@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize_ch, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
#Heatmap
pathways.show <- c("COLLAGEN")
par(mfrow=c(1,1))
netVisual_heatmap(cellchat_subset_ch, signaling = pathways.show, color.heatmap = "Reds")


# Chord diagram
#group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
#names(group.cellType) <- levels(cellchat_subset@idents)
#netVisual_chord_cell(cellchat_subset, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))
#> Plot the aggregated cell-cell communication network at the signaling pathway level
#contribution
cellchat_subset_ch <- netAnalysis_computeCentrality(cellchat_subset_ch)
netAnalysis_signalingRole_heatmap(cellchat_subset_ch, signaling = NULL, color.heatmap = "Reds", font.size = 5)

netAnalysis_contribution(cellchat_subset_ct, signaling = pathways.show)

#plotgeneExpression
plotGeneExpression(cellchat_subset_ch, signaling = "FN1", enriched.only = FALSE)

#contribution
netAnalysis_contribution(cellchat_subset_ch, signaling = pathways.show)

pairLR.CXCL <- extractEnrichedLR(cellchat_subset_ch, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair

# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat_subset_ch, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
#> [[1]]
# Circle plot
netVisual_individual(cellchat_subset_ch, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

###############
##regression###
##############
regression <- SetIdent(regression, value = regression@meta.data$new_clusters_id)
data.input_r <- as.data.frame(regression@assays$RNA@data)
cells_r <- regression@meta.data
cells_r <-cells_r[Cells(regression),]
data.input_r <- data.input_r[,row.names(cells_r)]
meta_r = regression@meta.data
cellchat_r <- createCellChat(object = as.matrix(data.input_r), meta = meta_r, group.by = "new_clusters_id")
cellchat_r@DB <- CellChatDB.use
cellchat_r <-subsetData(cellchat_r)
cellchat_r <- identifyOverExpressedGenes(cellchat_r)
cellchat_r <- identifyOverExpressedInteractions(cellchat_r)
cellchat_r <- computeCommunProb(cellchat_r)
cellchat_subset_r <- filterCommunication(cellchat_r, min.cells = 10)
cellchat_subset_r <- computeCommunProbPathway(cellchat_subset_r)
cellchat_subset_r <- aggregateNet(cellchat_subset_r)

groupSize_r <- as.numeric(table(cellchat_subset_r@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_subset_r@net$count, vertex.weight = groupSize_r, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_subset_r@net$weight, vertex.weight = groupSize_r, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

#Individuals
mat <- cellchat_subset_r@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize_r, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
#Heatmap
pathways.show <- c("NEGR")
par(mfrow=c(1,1))
netVisual_heatmap(cellchat_subset_r, signaling = pathways.show, color.heatmap = "Reds")


# Chord diagram
#group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
#names(group.cellType) <- levels(cellchat_subset@idents)
#netVisual_chord_cell(cellchat_subset, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))
#> Plot the aggregated cell-cell communication network at the signaling pathway level

#contribution
cellchat_subset_r <- netAnalysis_computeCentrality(cellchat_subset_r)
netAnalysis_signalingRole_heatmap(cellchat_subset_r, signaling = NULL, color.heatmap = "Reds", font.size = 6)
pathways.show <- c("NEGR")
netAnalysis_contribution(cellchat_subset_ct, signaling = pathways.show)

#plotgeneExpression
plotGeneExpression(cellchat_subset_ct, signaling = "FN1",  enriched.only = F)

#contribution
netAnalysis_contribution(cellchat_subset_ch, signaling = pathways.show)
pathways.show <- c( "COMPLEMENT")

pairLR.CXCL <- extractEnrichedLR(cellchat_subset_ch, signaling = pathways.show, geneLR.return = FALSE, enriched.only = F)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair

# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat_subset_ch, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
#> [[1]]
# Circle plot
netVisual_individual(cellchat_subset_ch, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")


#MERGE BACK
object.list_CT_CH <- list(CT = cellchat_subset_ct, CH = cellchat_subset_ch)
object.list_CH_R <- list(CH = cellchat_subset_ch, R= cellchat_subset_r)
object.list_CT_R <- list(CT = cellchat_subset_ct, R= cellchat_subset_r)
cellchat_CT_CH <- mergeCellChat(object.list_CT_CH, add.names = names(object.list_CT_CH))
cellchat_CH_R <- mergeCellChat(object.list_CH_R, add.names = names(object.list_CH_R))
cellchat_CT_R <- mergeCellChat(object.list_CT_R, add.names = names(object.list_CT_R))

p1 <-compareInteractions(cellchat_CT_CH, show.legend = F, group = c(1,2))
p2 <-compareInteractions(cellchat_CH_R, show.legend = F, group = c(1,2))
p3 <-compareInteractions(cellchat_CT_R, show.legend = F, group = c(1,2))
p1+p2+p3

gg1 <- netVisual_heatmap(cellchat_CT_CH)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat_CT_CH, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

#circle

weight.max_ct_ch <- getMaxWeight(object.list_CT_CH, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list_CT_CH)) {
  netVisual_circle(object.list_CT_CH[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max_ct_ch[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list_CT_CH)[i]))
}
cellchat_CT_CH <- netAnalysis_computeCentrality(cellchat_CT_CH)
netAnalysis_signalingRole_heatmap(cellchat_subset_r, signaling = NULL, color.heatmap = "Reds", font.size = 6)
pathways.show <- c("NEGR")
netAnalysis_contribution(cellchat_subset_ct, signaling = pathways.show)

#2D Space
num.link <- sapply(object.list_CT_CH, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list_CT_CH)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list_CT_CH[[i]], title = names(object.list_CT_CH)[i], weight.MinMax = weight.MinMax)
}

patchwork::wrap_plots(plots = gg)

#signalling associated to certain cell types

gg1 <- netAnalysis_signalingChanges_scatter(cellchat_CT_CH, idents.use = "HSCs")

gg2 <- netAnalysis_signalingChanges_scatter(cellchat_CT_R, idents.use = "HSCs")

patchwork::wrap_plots(plots = list(gg1,gg2))

#plotgeneExpression
plotGeneExpression(cellchat_CT_R, signaling = "ANGPTL", enriched.only = F, group.by = "orig.ident")
gg2 <-plotGeneExpression(cellchat_subset_r, signaling = "PTPRM",  enriched.only = F)

patchwork::wrap_plots(plots = list(gg1,gg2))

#plot gene expression in general
plotGeneExpression(cellchat_CT_R, signaling = "PTPRM",enriched.only = F , group.by = "orig.ident")

#Identify signaling groups based on their functional similarity
library(reticulate)
#py_install('umap-learn')
cellchat_CT_CH <- computeNetSimilarityPairwise(cellchat_CT_CH, type = "functional")
cellchat_CT_CH <- netEmbedding(cellchat_CT_CH, type = "functional")
cellchat_CT_CH <- netClustering(cellchat_CT_CH, type = "functional")
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat_CT_CH, type = "functional", label.size = 3.5)

cellchat_CH_R <- computeNetSimilarityPairwise(cellchat_CH_R, type = "functional")
cellchat_CH_R <- netEmbedding(cellchat_CH_R, type = "functional")
cellchat_CH_R <- netClustering(cellchat_CH_R, type = "functional")
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat_CH_R, type = "functional", label.size = 3.5)

cellchat_CT_R <- computeNetSimilarityPairwise(cellchat_CT_R, type = "functional")
cellchat_CT_R <- netEmbedding(cellchat_CT_R, type = "functional")
cellchat_CT_R <- netClustering(cellchat_CT_R, type = "functional")
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat_CT_R, type = "functional", label.size = 3.5)

#Identify signaling groups based on structure similarity
cellchat_CT_CH <- computeNetSimilarityPairwise(cellchat_CT_CH, type = "structural")
cellchat_CT_CH <- netEmbedding(cellchat_CT_CH, type = "structural")
cellchat_CT_CH <- netClustering(cellchat_CT_CH, type = "structural")
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat_CT_CH, type = "structural", label.size = 3.5)

cellchat_CH_R <- computeNetSimilarityPairwise(cellchat_CH_R, type = "structural")
cellchat_CH_R <- netEmbedding(cellchat_CH_R, type = "structural")
cellchat_CH_R <- netClustering(cellchat_CH_R, type = "structural")
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat_CH_R, type = "structural", label.size = 3.5)

cellchat_CT_R <- computeNetSimilarityPairwise(cellchat_CT_R, type = "structural")
cellchat_CT_R <- netEmbedding(cellchat_CT_R, type = "structural")
cellchat_CT_R <- netClustering(cellchat_CT_R, type = "structural")
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat_CH_R, type = "structural", label.size = 3.5)

#Compute and visualize the pathway distance in the learned joint manifold
rankSimilarity(cellchat_CT_CH, type = "functional")
rankSimilarity(cellchat_CH_R, type = "functional")
rankSimilarity(cellchat_CT_R, type = "functional")

#Compare the overall information flow of each signaling pathway
rankNet(cellchat_CT_CH, mode = "comparison", stacked = T, do.stat = TRUE, font.size = 6)
rankNet(cellchat_CH_R, mode = "comparison", stacked = T, do.stat = TRUE, font.size = 6)
rankNet(cellchat_CT_R, mode = "comparison", stacked = T, do.stat = TRUE, font.size = 6)

#Identify dysfunctional signaling by comparing the communication probabities
netVisual_bubble(cellchat_CH_R, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), angle.x = 45, font.size = 6)
gg1 <- netVisual_bubble(cellchat_CT_R, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in R", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat_CT_R, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in R", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2

#Identify dysfunctional signaling by using differential expression analysis
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "R"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat_CT_R <- identifyOverExpressedGenes(cellchat_CT_R, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat_CT_R, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat_CT_R, net = net, datasets = "R",ligand.logFC = 0.2, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat_CT_R, net = net, datasets = "CT",ligand.logFC = -0.1, receptor.logFC = -0.1)

#Individual genens
gene.up <- extractGeneSubsetFromPair(net.up, cellchat_CT_R)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat_CT_R)

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat_CT_R, pairLR.use = pairLR.use.up, sources.use = 4, targets.use = c(5:11), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list_CT_R)[2]))

#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat_CT_R, pairLR.use = pairLR.use.down, sources.use = 4, targets.use = c(5:11), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list_CT_R)[2]))
#> Comparing communications on a merged object
gg1 + gg2

#Visualize the upgulated and down-regulated signaling ligand-receptor pairs using Chord diagram
# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list_CT_R[[2]], sources.use = 4, targets.use = c(5:11), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list_CT_R)[2]))
netVisual_chord_gene(object.list_CT_R[[1]], sources.use = 4, targets.use = c(5:11), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list_CT_R)[2]))

# visualize the enriched ligands in the first condition
computeEnrichmentScore(net.down, species = 'human') #no funciona

#Visually compare cell-cell communication using Hierarchy plot, Circle plot or Chord diagram

pathways.show <- c("APP") 
weight.max <- getMaxWeight(object.list_CT_CH, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list_CT_CH)) {
  netVisual_aggregate(object.list_CT_CH[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, thresh = 0.9, signaling.name = paste(pathways.show, names(object.list_CT_CH)[i]),)
}

ht <- list()
for (i in 1:length(object.list_CT_CH)) {
  ht[[i]] <- netVisual_heatmap(object.list_CT_CH[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list_CT_CH)[i]))
}

#> Do heatmap based on a single object
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

#Compare outgoing (or incoming) signaling associated with each cell population

i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list_CT_R[[i]]@netP$pathways, object.list_CT_R[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list_CT_R[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list_CT_R)[i], width = 5, height = 15, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list_CT_R[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list_CT_R)[i+1], width = 5, height = 15, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

#Compare the signaling gene expression distribution between different datasets

cellchat_CT_CH@meta$datasets = factor(cellchat_CT_CH@meta$datasets, levels = c("CT", "CH")) # set factor level
plotGeneExpression(cellchat_CT_CH, signaling = "CADM", split.by = "datasets", colors.ggplot = T, enriched.only = F)

##IMMUNE SEURAT OBJECT##

data.input_immune <- as.data.frame(immune@assays$RNA@data)
cells <- all.liver.combined@meta.data
cells <-cells[Cells(all.liver.combined),]
data.input <- data.input[,row.names(cells)]
meta = all.liver.combined@meta.data
cellchat <- createCellChat(object = as.matrix(data.input), meta = meta, group.by = "new_clusters_id")
cellchat@DB <- CellChatDB.use
cellchat <-subsetData(cellchat)
future::plan("multiprocess", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#optional cellchat2 <- projectData(cellchat2, PPI.human)
cellchat <- computeCommunProb(cellchat)
cellchat_subset <- filterCommunication(cellchat, min.cells = 10)
cellchat_subset <- computeCommunProbPathway(cellchat_subset)
cellchat_subset <- aggregateNet(cellchat_subset)

