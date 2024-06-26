# Librerias
library(dplyr)
library(Seurat)
library(tidyverse)
library(gridExtra)

# Datos procesados
seurat.integrated <- readRDS("~/sc-RNAseq1/analisis/1calidadQC/seurat_integrated.rds")
DimPlot(seurat.integrated, reduction = "umap", group.by = c("Condition", "seurat_clusters"),label = T)

# Marcadores monocitos y celulas dendriticas
genes_mono<-c("CTSS", "FCN1", "NEAT1", "LYZ", "PSAP", "S100A9", "AIF1", "MNDA", "SERPINA1", "TYROBP")
genes_DC<-c("CD74", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "CCDC88A", "HLA-DRA", "HLA-DMA", "CST3", "HLA-DQB1", "HLA-DRB1")

FeaturePlot(seurat.integrated, reduction = "umap", features = genes_mono)
FeaturePlot(seurat.integrated, reduction = "umap", features = genes_DC)

# Marcadores linfocitos T y NK
genes_CD4positive<-c("IL7R", "MAL", "LTB", "CD4", "LDHB", "TPT1", "TRAC", "TMSB10", "CD3D", "CD3G")
genes_CD8positive<-c("CD8B", "CD8A", "CD3D", "TMSB10", "HCST", "CD3G", "LINC02446", "CTSW", "CD3E", "TRAC")
genes_otherT<-c("CD3D", "TRDC", "GZMK", "KLRB1", "NKG7", "TRGC2", "CST7", "LYAR", "KLRG1", "GZMA")
genes_NK<-c("NKG7", "KLRD1", "TYROBP", "GNLY", "FCER1G", "PRF1", "CD247", "KLRF1", "CST7", "GZMB")

FeaturePlot(seurat.integrated, reduction = "umap", features = genes_CD4positive)
FeaturePlot(seurat.integrated, reduction = "umap", features = genes_CD8positive)
FeaturePlot(seurat.integrated, reduction = "umap", features = genes_otherT)
FeaturePlot(seurat.integrated, reduction = "umap", features = genes_NK)

# Marcadores linfocitos B
genes_B<-c("CD79A", "RALGPS2","CD79B", "MS4A1", "BANK1", "CD74", "TNFRSF13C", "HLA-DQA1", "IGHM", "MEF2C")

FeaturePlot(seurat.integrated, reduction = "umap", features = genes_B)

# Marcadores plaquetas y plasmablastos
genes_platelet <- c("PPBP", "PF4", "NRGN", "GNG11", "CAVIN2", "TUBB1", "CLU", "HIST1H2AC", "RGS18", "GP9")
genes_plasmablast<-c("IGHA2", "MZB1", "TNFRSF17", "DERL3", "TXNDC5", "TNFRSF13B", "POU2AF1", "CPNE5", "NT5DC2")

plot.list<-list() 
for (i in c(genes_platelet)){
  plot.list[[i]]<-DotPlot(seurat.integrated,group.by = "RNA_snn_res.0.5",features = c(i))+ theme(legend.position = 'none')
  
}
do.call(grid.arrange, c(plot.list[1:10],ncol=10))

plot.list<-list() 
for (i in c(genes_plasmablast)){
  plot.list[[i]]<-DotPlot(seurat.integrated,group.by = "RNA_snn_res.0.5",features = c(i))+ theme(legend.position = 'none')
  
}
do.call(grid.arrange, c(plot.list[1:9],ncol=9))
