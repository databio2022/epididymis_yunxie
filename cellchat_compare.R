#!/usr/bin/Rscript

library(CellChat)
library(patchwork)
library(ggplot2)
library(igraph)  

data.dir <- 'cellchat_diff' 
dir.create(data.dir)

input_filepath=read.table("experiment.txt",header=TRUE,sep="\t")# experiment config file

object.list=list()
for(i in 1:nrow(input_filepath)){
  object.list[[input_filepath$group[i]]]=readRDS(input_filepath$rds[i])
}

cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat

setwd(data.dir)

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
pdf(file="compare_interactions_circleplot.pdf",width=16,height=8)
par(mfrow = c(1,length(object.list)), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

weight.max <- getMaxWeight(object.list, attribute = c("idents","weight")) 
pdf(file="compare_interactions_strength_circleplot.pdf",width=16,height=8)
par(mfrow = c(1,length(object.list)), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$weight, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Interaction weights/strength - ", names(object.list)[i]))
}
dev.off()


#outgoing and incoming
library(ComplexHeatmap)
pathway.union=c()
for(i in 1:length(object.list)){
  pathway.union <- append(pathway.union,object.list[[i]]@netP$pathways)
}

pattern = c("outgoing", "incoming")
for(k in pattern){
p = netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = k, signaling = pathway.union, title = names(object.list)[1], width = 10, height = 20)#第一个对象热图

for(i in 2:length(object.list)){
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = k, signaling = pathway.union, title = names(object.list)[i], width = 10, height = 20)#第一个对象热图
p=p+ht1
}

pdf(file=paste0("compare_",k,".pdf"),width=25,height=25)
draw(p, ht_gap = unit(0.5, "cm"))
dev.off()
}


saveRDS(cellchat, file = "cellchat_diff_analysis.rds")
