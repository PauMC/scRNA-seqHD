# Librerias
library(dplyr)
library(Seurat)
library(tidyverse)
library(gridExtra)

# Datos procesados
seurat.integrated <- readRDS("~/sc-RNAseq1/analisis/1calidadQC/seurat_integrated.rds")
DimPlot(seurat.integrated, reduction = "umap", group.by = c("Condition", "seurat_clusters"),label = T)

# Establecer Ident
seurat.integrated<-SetIdent(seurat.integrated,value = "seurat_clusters")

# Seleccionamos los clusters de interes
seuratB<-subset(seurat.integrated,idents = c("10","11","12","16"))

seuratB[["RNA"]] <- split(seuratB[["RNA"]], f = seuratB$Condition)

seuratB <- FindVariableFeatures(seuratB)
seuratB <- ScaleData(seuratB)
seuratB <- RunPCA(seuratB)
ElbowPlot(seuratB)

seuratB <- IntegrateLayers(object = seuratB, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                           verbose = FALSE)

seuratB[["RNA"]] <- JoinLayers(seuratB[["RNA"]])
seuratB <- FindNeighbors(seuratB, reduction = "integrated.cca", dims = 1:20)
seuratB <- FindClusters(seuratB, resolution = 0.4)
seuratB <- RunUMAP(seuratB, dims = 1:20, reduction = "integrated.cca")

DimPlot(seuratB, reduction = "umap", group.by = c("Condition", "seurat_clusters"),label = T)

# Marcadores linfocitos B
genes_Bintermediate<-c("TNFRSF13B","LINC01857")
genes_Bmemory<-c("COCH","SSPN","TEX9","TNFRSF13C","LINC01781")
genes_Bnaive<-c("IL4R","CXCR4","BTG1","TCL1A","YBX3")

plot.list<-list() 
for (i in c(genes_Bintermediate,genes_Bmemory,genes_Bnaive)){
  plot.list[[i]]<-DotPlot(seuratB,group.by = "RNA_snn_res.0.4",features = c(i))+ theme(legend.position = 'none')
  
}
do.call(grid.arrange, c(plot.list[1:12],ncol=12))

# Marcadores plasmablastos
genes_plasmablast<-c("IGHA2", "MZB1", "TNFRSF17", "DERL3", "TXNDC5", "TNFRSF13B", "POU2AF1", "CPNE5", "NT5DC2")

plot.list<-list() 
for (i in c(genes_plasmablast)){
  plot.list[[i]]<-DotPlot(seuratB,group.by = "RNA_snn_res.0.4",features = c(i))+ theme(legend.position = 'none')
  
}
do.call(grid.arrange, c(plot.list[1:9],ncol=9))

# Marcadores MAIT
genes_MAIT<-c("KLRB1", "NKG7", "GZMK", "IL7R", "SLC4A10", "GZMA", "CXCR6", "PRSS35", "RBM24", "NCR3")
plot.list<-list() 
for (i in c(genes_MAIT)){
  plot.list[[i]]<-DotPlot(seuratB,group.by = "RNA_snn_res.0.4",features = c(i))+ theme(legend.position = 'none')
  
}
do.call(grid.arrange, c(plot.list[1:10],ncol=10))

# Porcentaje mitocondrial
VlnPlot(seuratB, features = "percent.mt")

# Renombrar clusters
new.cluster.ids <- c("Naive", "Memory", "Memory", "Naive", "Intermediate", "Doblete B-plasmablasto", "Doblete B-MAIT")

names(new.cluster.ids) <- levels(seuratB)

seuratB <- RenameIdents(seuratB, new.cluster.ids)

plot <- DimPlot(seuratB, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
plot + coord_cartesian(xlim = c(-15, 6))

saveRDS(seuratB,"seuratB.rds")

