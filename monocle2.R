#pseudotime analysis
library(monocle)
library(Seurat)
library(RColorBrewer)

seurat=readRDS(file= "seuratObject.rds")
rootcelltype="Basal cells"
matrix=seurat@assays$RNA@counts

matrix=matrix[which(rowSums(matrix) > 10),]
write.csv(matrix,file='matrix.csv')
#featuredata
feature_ann<-data.frame(gene_id=rownames(matrix),gene_short_name=rownames(matrix))
rownames(feature_ann)<-rownames(matrix)
fd<-new("AnnotatedDataFrame", data = feature_ann)
#metadata to phenodata
sample_ann<-seurat@meta.data
rownames(sample_ann)<-colnames(matrix)
pd<-new("AnnotatedDataFrame", data =sample_ann)
cds<-newCellDataSet(matrix,phenoData =pd,featureData =fd,expressionFamily=negbinomial.size())

head(pData(cds))
head(fData(cds))
cds <- estimateSizeFactors(cds) 
cds <- estimateDispersions(cds) 
cds <- detectGenes(cds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(cds),num_cells_expressed >= 10))
diff_celltype <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr = "~celltype",cores=12)
diff_celltype<- diff_celltype[order(diff_celltype$qval),]
write.csv(diff_celltype, file = "diff_celltype.csv")
ordering_genes <- row.names(diff_celltype[1:1000,])
cds <- setOrderingFilter(cds,ordering_genes = ordering_genes)
cds <- reduceDimension(cds, method = 'DDRTree')
cds <- orderCells(cds)
GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$celltype)[,rootcelltype]
    return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
    } else {return (1)}}
cds <- orderCells(cds, root_state = GM_state(cds))


pdf("celltype_traj.pdf")
plot_cell_trajectory(cds, color_by = "celltype")
dev.off()

pdf("celltype_split_traj.pdf")
plot_cell_trajectory(cds, color_by = "celltype") +
  facet_wrap(~celltype, nrow = 2)
dev.off()

###pseudotime dependent gene
diff_test_res=differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr = "~sm.ns(Pseudotime)",cores=12)

##top30
top30g=rownames(diff_test_res[order(diff_test_res$qval),])[1:30]
pdf("top30_genes_in_pseudotime.pdf")
plot_genes_in_pseudotime(cds[top30g,], color_by="Pseudotime", ncol = 3)+scale_color_viridis_c()
dev.off()

saveRDS(cds, file = "monocle2.rds")
