# Librerias
library(dplyr)
library(Seurat)
library(tidyverse)
library(gridExtra)

# Datos procesados
seurat.integrated <- readRDS("~/sc-RNAseq1/analisis/1calidadQC/seurat_integrated.rds")

# Establecer Ident
seurat.integrated<-SetIdent(seurat.integrated,value = "seurat_clusters")

# Seleccionamos los clusters de interés
seuratTyNK<-subset(seurat.integrated,idents = c("1","2","3","4","5","6","7","9"))

seuratTyNK[["RNA"]] <- split(seuratTyNK[["RNA"]], f = seuratTyNK$Condition)

seuratTyNK <- FindVariableFeatures(seuratTyNK)
seuratTyNK <- ScaleData(seuratTyNK)
seuratTyNK <- RunPCA(seuratTyNK)
ElbowPlot(seuratTyNK)

seuratTyNK <- IntegrateLayers(object = seuratTyNK, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                              verbose = FALSE)

seuratTyNK[["RNA"]] <- JoinLayers(seuratTyNK[["RNA"]])
seuratTyNK <- FindNeighbors(seuratTyNK, reduction = "integrated.cca", dims = 1:20)
seuratTyNK <- FindClusters(seuratTyNK, resolution = 0.5)
seuratTyNK <- RunUMAP(seuratTyNK, dims = 1:20, reduction = "integrated.cca")
DimPlot(seuratTyNK, reduction = "umap", group.by = c("Condition", "seurat_clusters"),label = T)

# Marcadores linfocitos T
genes_linfT<- c("CD3D","CD3E","CD3G")

library(gridExtra)
plot.list<-list() 
for (i in c(genes_linfT)){
  plot.list[[i]]<-DotPlot(seuratTyNK,group.by = "RNA_snn_res.0.5",features = c(i))+ theme(legend.position = 'none')
  
}
do.call(grid.arrange, c(plot.list[1:3],ncol=3))

genes_CD4positive<-c("IL7R", "MAL", "LTB", "CD4", "LDHB", "TPT1", "TMSB10")
genes_CD8positive<-c("CD8B", "CD8A", "TMSB10", "HCST","LINC02446", "CTSW")
plot.list<-list() 
for (i in c(genes_CD4positive,genes_CD8positive)){
  plot.list[[i]]<-DotPlot(seuratTyNK,group.by = "RNA_snn_res.0.5",features = c(i))+ theme(legend.position = 'none')
  
}
do.call(grid.arrange, c(plot.list[1:13],ncol=13))

genes_CD8naive<-c("CD8B", "S100B", "CCR7", "RGS10", "NOSIP", "LINC02446", "LEF1", "CRTAM", "CD8A", "OXNAD1")
plot.list<-list() 
for (i in c(genes_CD8naive)){
  plot.list[[i]]<-DotPlot(seuratTyNK,group.by = "RNA_snn_res.0.5",features = c(i))+ theme(legend.position = 'none')
  
}
do.call(grid.arrange, c(plot.list[1:10],ncol=10))

# Marcadores NK
genes_NK<-c("NKG7", "KLRD1", "TYROBP", "GNLY", "FCER1G", "PRF1", "CD247", "KLRF1", "CST7", "GZMB")

plot.list<-list() 
for (i in c(genes_NK)){
  plot.list[[i]]<-DotPlot(seuratTyNK,group.by = "RNA_snn_res.0.5",features = c(i))+ theme(legend.position = 'none')
  
}
do.call(grid.arrange, c(plot.list[1:10],ncol=10))

# Marcadores MAIT
genes_MAIT<-c("KLRB1", "NKG7", "GZMK", "IL7R", "SLC4A10", "GZMA", "CXCR6", "PRSS35", "RBM24", "NCR3")
plot.list<-list() 
for (i in c(genes_MAIT)){
  plot.list[[i]]<-DotPlot(seuratTyNK,group.by = "RNA_snn_res.0.5",features = c(i))+ theme(legend.position = 'none')
  
}
do.call(grid.arrange, c(plot.list[1:10],ncol=10))

# Porcentaje mitocondrial
VlnPlot(seuratTyNK, features = "percent.mt")

# https://satijalab.org/seurat/articles/pbmc3k_tutorial 
# Renombrar clusters
new.cluster.ids <- c("CD4", "CD4", "CD8", "Baja calidad", "CD8", " Baja calidad ", "CD8 naive", "NK", "CD4", "NK", "MAIT", "NK")

names(new.cluster.ids) <- levels(seuratTyNK)

seuratTyNK <- RenameIdents(seuratTyNK, new.cluster.ids)

DimPlot(seuratTyNK, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

saveRDS(seuratTyNK,"seuratTyNK.rds")
