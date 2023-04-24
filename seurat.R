#/user/bin/Rscript
#integrate analysis for mulitple single cell sequencing samples

dir.create("outFigures")
dir.create("outData")

library(Seurat)
library(cowplot)
library(patchwork)
library(Matrix)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggthemes)

experi=read.table("experiment.txt",header=TRUE,sep="\t",check.names=FALSE) #experiment config file
allsample=NULL

for(i in 1:length(experi$filePath)){
	d=read.table(experi$filePath[i],header=TRUE,sep="\t")
	d=d[!duplicated(d[,1]),]
	rownames(d)=d[,1]
	d=d[,-1]
	d=as.matrix(d)
	sparse <- Matrix(d, sparse = T )
	data=sparse
	# Initialize the Seurat object with the raw (non-normalized data).
	data <- CreateSeuratObject(counts = data, min.cells = 3, min.features = 200)
	# The [[ operator can add columns to object metadata. This is a great place to stash QC stats

	data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^mt-") #mouse
	
	for(j in 2:length(names(experi))){
		data[[names(experi)[j]]]=experi[i,j]
	}
	# Visualize QC metrics as a violin plot
	QC1=VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
	ggsave(paste0("outFigures/",experi$sample[i],".QC1.pdf"),QC1,width=16,height=9)
	# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
	# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
	plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
	plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
	QC2=CombinePlots(plots = list(plot1, plot2))
	ggsave(paste0("outFigures/",experi$sample[i],".QC2.pdf"),QC2,width=16,height=9)
	#filter cells
	data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 20)
	allsample=append(allsample,data)
}

alls1=allsample
alls1[[1]]=NULL
alls=merge(x=allsample[[1]],y=alls1)

ifnb.list <- SplitObject(alls, split.by = "sample")

ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = ifnb.list)

anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)

combined <- IntegrateData(anchorset = anchors)


DefaultAssay(combined) <- "integrated"

# Run the standard workflow for visualization and clustering
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
# UMAP and Clustering
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.5)

p1 <- DimPlot(combined, reduction = "umap", group.by = "group")
p2 <- DimPlot(combined, reduction = "umap", label = TRUE, repel = TRUE)
pdf(file="outFigures/seurat_integrated_umap.pdf",width=16)
plot_grid(p1, p2)
dev.off()

pdf(file="outFigures/seurat_integrated_umap_splitgroup.pdf",width=16)
print(DimPlot(combined, reduction = "umap", split.by = "group"))
dev.off()


DefaultAssay(combined) <- "RNA"


marker_table=NULL

all_diff_gene=NULL

# find markers for every cluster compared to all remaining cells, report only the positive
markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers=markers[markers$p_val_adj <= 0.05,]
all_diff_gene=markers
marker_table=markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) # top 10


write.table(marker_table,"outData/top10marker.txt",quote = FALSE,sep="\t",col.names=NA)
write.table(all_diff_gene,"outData/all.diffgene_onlypos_padj005.txt",quote = FALSE,sep="\t",col.names=NA)


# cell annotation and differential expression analysis
cellanno=read.table("cell-annotation.txt",sep="\t",header=TRUE) #manual annotation results
cell_label=NULL
for(i in 1:length(levels(combined))){
	cells.use <- WhichCells(combined, idents = as.character(cellanno$clusters)[i])
	combined <- SetIdent(combined, cells = cells.use, value = as.character(cellanno$labels)[i])
}

p1 <- DimPlot(combined, reduction = "umap", group.by = "group")
p2 <- DimPlot(combined, reduction = "umap", label = TRUE, repel = TRUE)
pdf(file="outFigures/seurat_integrated_cellann.pdf",width=16)
plot_grid(p1, p2)
dev.off()

pdf(file="outFigures/seurat_integrated_cellann_splitgroup.pdf",width=16)
print(DimPlot(combined, reduction = "umap", split.by = "group",label = TRUE, repel = TRUE))
dev.off()

#Identify differential expressed genes for same cell type across conditions
combined$celltype.group <- paste(Idents(combined), combined$group, sep = "_")
combined$celltype <- Idents(combined)

Idents(combined) <- "celltype.group"

experi=read.table("vs.txt",header=TRUE,sep="\t",check.names=FALSE) # experiment comparsion settings file
for(i in unique(cellanno$labels)){
	for(ca in experi$vs){
		if(ca != "no"){
			cas=strsplit(ca,"_")    
			if(length(which(combined@meta.data$group==cas[[1]][1] & combined@meta.data$celltype==i)) > 2 && length(which(combined@meta.data$group==cas[[1]][3] & combined@meta.data$celltype==i)) > 2){
				diff <- FindMarkers(combined, ident.1 = paste0(i,"_",cas[[1]][1]), ident.2 = paste0(i,"_",cas[[1]][3]), verbose = FALSE)
				write.table(diff,paste0("outData/",i,"_",ca,"_deanalysis.xls"),quote = FALSE,sep="\t",col.names=NA)
			}
		}
	}
}

####heatmap####

markers=all_diff_gene  # only pos padj005

###2fc###
combined <- ScaleData(combined, features = unique(markers[markers$avg_log2FC >= 1,]$gene)) #2fc
Idents(combined)=combined$celltype
####average heatmap####
cluster.averages <- AverageExpression(combined, return.seurat = TRUE)
p=DoHeatmap(object = cluster.averages,features=markers[markers$avg_log2FC >= 1,]$gene,size = 3,draw.lines = FALSE)+theme(axis.text.y = element_blank())

ggsave("outFigures/celltype_2fc_markers_heatmap_averageexp.pdf",p,width=8,height=12)

saveRDS(combined,"outData/seuratObject.rds")
