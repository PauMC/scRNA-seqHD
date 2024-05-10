
#Librerias
library(dplyr)
library(Seurat)
library(tidyverse)
library(gridExtra)
library(ggrepel)

# Función para encontrar marcadores y escribir los resultados en un archivo CSV
find_and_write_markers <- function(subset_seurat, idents, file_name) {
  Idents(subset_seurat) <- "Condition"
  de_markers <- FindMarkers(subset_seurat, ident.1 = "hd", ident.2 = "control", slot = "counts", test.use = "MAST", verbose = FALSE, latent.vars = "sex")
  de_markers$gene <- rownames(de_markers)
  write.csv(de_markers, file_name)
}

# LinfB
seuratB <- readRDS("~/sc-RNAseq1/analisis/2clustering/seuratB.rds")

DimPlot(seuratB, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

countscells<-data.frame(seuratB$RNA$counts)
filter<-Features(seuratB)[Features(seuratB) %in% c("XIST")]
sexmarkers<-data.frame(t(countscells[rownames(countscells) %in% filter,]))

seuratB$sex <- ifelse(sexmarkers$XIST>0,"Mujer","Hombre")

identities <- c("Naive", "Memory", "Intermediate")
file_names <- c("de_markers_seuratB_naive_MAST.csv", "de_markers_seuratB_memory_MAST.csv", "de_markers_seuratB_intermediate_MAST.csv")

for (i in seq_along(identities)) {
  subset_seurat <- subset(x = seuratB, idents = identities[i])
  find_and_write_markers(subset_seurat, identities[i], file_names[i])
}

# MonoDC
seuratmonoDC <- readRDS("~/sc-RNAseq1/analisis/2clustering/seuratmonoDC.rds")

DimPlot(seuratmonoDC, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

countscells<-data.frame(seuratmonoDC$RNA$counts)
filter<-Features(seuratmonoDC)[Features(seuratmonoDC) %in% c("XIST")]
sexmarkers<-data.frame(t(countscells[rownames(countscells) %in% filter,]))

seuratmonoDC$sex <- ifelse(sexmarkers$XIST>0,"Mujer","Hombre")

identities <- c("CD14", "CD16", "cDC2")
file_names <- c("de_markers_monoDC_CD14_MAST.csv", "de_markers_monoDC_CD16_MAST.csv", "de_markers_monoDC_cDC2_MAST.csv")

for (i in seq_along(identities)) {
  subset_seurat <- subset(x = seuratmonoDC, idents = identities[i])
  find_and_write_markers(subset_seurat, identities[i], file_names[i])
}

# LinfTyNK
seuratTyNK <- readRDS("~/sc-RNAseq1/analisis/2clustering/seuratTyNK.rds")

DimPlot(seuratTyNK, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

countscells<-data.frame(seuratTyNK$RNA$counts)
filter<-Features(seuratTyNK)[Features(seuratTyNK) %in% c("XIST")]
sexmarkers<-data.frame(t(countscells[rownames(countscells) %in% filter,]))

seuratTyNK$sex <- ifelse(sexmarkers$XIST>0,"Mujer","Hombre")

identities <- c("CD4", "CD8", "NK", "CD8 naive", "MAIT")
file_names <- c("de_markers_linfTyNK_CD4_MAST.csv", "de_markers_linfTyNK_CD8_MAST.csv", "de_markers_linfTyNK_NK_MAST.csv",
                     "de_markers_linfTyNK_CD8_naive_MAST.csv", "de_markers_linfTyNK_MAIT_MAST.csv")

for (i in seq_along(identities)) {
  subset_seurat <- subset(x = seuratTyNK, idents = identities[i])
  find_and_write_markers(subset_seurat, identities[i], file_names[i])
}

