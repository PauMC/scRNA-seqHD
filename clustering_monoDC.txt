# Librerias
library(dplyr)
library(Seurat)
library(tidyverse)
library(gridExtra)

# Datos procesados
seurat.integrated <- readRDS("~/sc-RNAseq1/analisis/1calidadQC/seurat_integrated.rds")

# Establecer Ident
seurat.integrated<-SetIdent(seurat.integrated,value = "seurat_clusters")

# Seleccionamos los clusters de interes
seuratmonoDC<-subset(seurat.integrated,idents = c("0","8","13","14"))

seuratmonoDC[["RNA"]] <- split(seuratmonoDC[["RNA"]], f = seuratmonoDC$Condition)

seuratmonoDC <- FindVariableFeatures(seuratmonoDC)
seuratmonoDC <- ScaleData(seuratmonoDC)
seuratmonoDC <- RunPCA(seuratmonoDC)
ElbowPlot(seuratmonoDC)

seuratmonoDC <- IntegrateLayers(object = seuratmonoDC, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                                verbose = FALSE)

seuratmonoDC[["RNA"]] <- JoinLayers(seuratmonoDC[["RNA"]])
seuratmonoDC <- FindNeighbors(seuratmonoDC, reduction = "integrated.cca", dims = 1:20)
seuratmonoDC <- FindClusters(seuratmonoDC, resolution = 0.4)
seuratmonoDC <- RunUMAP(seuratmonoDC, dims = 1:20, reduction = "integrated.cca")
DimPlot(seuratmonoDC, reduction = "umap", group.by = c("Condition", "seurat_clusters"),label = T)

# Marcadores monocitos CD14 y CD16
genes_CD14Mono<-c("S100A9", "CTSS", "S100A8", "LYZ", "VCAN", "S100A12", "IL1B", "CD14", "G0S2", "FCN1")
genes_CD16Mono<-c("CDKN1C", "FCGR3A", "PTPRC", "LST1", "IER5", "MS4A7", "RHOC", "IFITM3", "AIF1", "HES4")

library(gridExtra)
plot.list<-list() 
for (i in c(genes_CD14Mono,genes_CD16Mono)){
  plot.list[[i]]<-DotPlot(seuratmonoDC,group.by = "RNA_snn_res.0.4",features = c(i))+ theme(legend.position = 'none')
  
}
do.call(grid.arrange, c(plot.list[1:21],ncol=21))

# Marcadores celulas MAIT
genes_MAIT<-c("KLRB1", "NKG7", "GZMK", "IL7R", "SLC4A10", "GZMA", "CXCR6", "PRSS35", "RBM24", "NCR3")
plot.list<-list() 
for (i in c(genes_MAIT)){
  plot.list[[i]]<-DotPlot(seuratmonoDC,group.by = "RNA_snn_res.0.4",features = c(i))+ theme(legend.position = 'none')
  
}
do.call(grid.arrange, c(plot.list[1:9],ncol=9))

# Marcadores plasmablastos
genes_plasmablast<-c("IGHA2", "MZB1", "TNFRSF17", "DERL3", "TXNDC5", "TNFRSF13B", "POU2AF1", "CPNE5", "NT5DC2")

plot.list<-list() 
for (i in c(genes_plasmablast)){
  plot.list[[i]]<-DotPlot(seuratmonoDC,group.by = "RNA_snn_res.0.4",features = c(i))+ theme(legend.position = 'none')
  
}
do.call(grid.arrange, c(plot.list[1:9],ncol=9))

# Marcadores celulas dendriticas cDC2
genes_cDC2<-c("FCER1A", "HLA-DQA1", "CLEC10A", "CD1C", "ENHO", "PLD4", "GSN", "SLC38A1", "NDRG2", "AFF3")
plot.list<-list() 
for (i in c(genes_cDC2)){
  plot.list[[i]]<-DotPlot(seuratmonoDC,group.by = "RNA_snn_res.0.4",features = c(i))+ theme(legend.position = 'none')
  
}
do.call(grid.arrange, c(plot.list[1:10],ncol=10))

# Porcentaje mitocondrial
VlnPlot(seuratmonoDC, features = "percent.mt")

# https://satijalab.org/seurat/articles/pbmc3k_tutorial 
# Renombramos clusters
new.cluster.ids <- c("CD14", "CD14", "CD14", "CD14", "CD16", "CD16", "cDC2", "Doblete mieloide-MAIT", "Doblete mieloide-plasmablasto")

names(new.cluster.ids) <- levels(seuratmonoDC)

seuratmonoDC <- RenameIdents(seuratmonoDC, new.cluster.ids)

plot <- DimPlot(seuratmonoDC, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
plot + coord_cartesian(xlim = c(-25, 6))

saveRDS(seuratmonoDC,"seuratmonoDC.rds")

