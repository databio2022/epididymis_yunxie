#!/user/bin/Rscript
# for single dataset（young or old）

library(Seurat)
library(ggplot2)
library(CellChat)
library(patchwork)
library(NMF)


options(stringsAsFactors = FALSE)

seurat=readRDS(file= "seuratObject.rds")

data.input = as.matrix(GetAssayData(seurat, assay = "RNA", slot = "data"))
meta = seurat@meta.data # a dataframe with rownames containing cell mata data
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")

levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

CellChatDB <- CellChatDB.mouse

CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- projectData(cellchat, PPI.mouse)


#################################Inference of cell-cell communication network#####################################
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat@netP$pathways

saveRDS(cellchat, file = "cellchat_single.rds")
