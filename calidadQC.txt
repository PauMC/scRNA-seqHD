# Librerias necesarias
library(dplyr)
library(Seurat)
library(tidyverse)
library(scDblFinder)

# Objecto Seurat
## Cargamos nuestras muestras de PBMC
control.data<-Read10X(data.dir = "/datos/paula/sc-RNAseq1/run_cellranger_count_control/outs/filtered_feature_bc_matrix/")
hd.data<-Read10X(data.dir = "/datos/paula/sc-RNAseq1/run_cellranger_count_hd/outs/filtered_feature_bc_matrix/")

## Creamos el objeto Seurat
control<-CreateSeuratObject(counts = control.data, min.cells = 3, min.features = 200)
hd<-CreateSeuratObject(counts = hd.data, min.cells = 3, min.features = 200)

## Juntamos los datos 
merged_seurat <- merge(control, y = hd,
                       add.cell.ids = c("control","hd"),
                       project = 'HD', merge.data = TRUE)

Layers(merged_seurat[["RNA"]])
merged_seurat[["RNA"]] <- JoinLayers(merged_seurat[["RNA"]])

View(merged_seurat@meta.data)

## Creamos una columna de muestra en el metadata
merged_seurat$sample <- rownames(merged_seurat@meta.data)

## Separamos la columna
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample', into = c('Condition', 'Barcode'), sep = '_')

# QC 
## Añadimos una columna al metadata con el porcentaje mitocondrial
merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "^MT-")

## Visualización de las métricas del QC
VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "Condition")

## Visualización de las relaciones del número de moléculas con el porcentaje mitocondrial o el número de genes
plot1 <- FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Condition")
plot2 <- FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Condition")
plot1 + plot2

# Filtrado
merged_seurat <- subset(merged_seurat, subset = nFeature_RNA > 800 & nFeature_RNA < 6500 & percent.mt < 20)
merged_seurat

# Normalizado
merged_seurat <- NormalizeData(merged_seurat)

# Cálculo del número de dobletes
dbl.dens <- computeDoubletDensity(as(merged_seurat[["RNA"]]$data, Class = "dgCMatrix"))
summary(dbl.dens)
merged_seurat$doublet_score <- dbl.dens
plot(merged_seurat$doublet_score)
merged_seurat$doublet_dblfinder <- merged_seurat$doublet_score > 5
table(merged_seurat$doublet_dblfinder)
merged_seurat <- subset(merged_seurat, doublet_score < 5)

# Análisis sin integrar los datos
merged_seurat <- FindVariableFeatures(merged_seurat)
merged_seurat <- ScaleData(merged_seurat)
merged_seurat <- RunPCA(merged_seurat)
ElbowPlot(merged_seurat)

# Reducción de dimensionalidad lineal PCA
merged_seurat <- FindNeighbors(merged_seurat, dims = 1:20, reduction = "pca")
merged_seurat <- FindClusters(merged_seurat, resolution = 0.5, cluster.name = "unintegrated_clusters")
merged_seurat <- RunUMAP(object = merged_seurat, dims = 1:20, reduction.name = "umap.unintegrated")

DimPlot(merged_seurat, reduction = "umap.unintegrated", group.by = c("Condition", "seurat_clusters"),label = T)

# Integrar los datos
merged_seurat[["RNA"]] <- split(merged_seurat[["RNA"]], f = merged_seurat$Condition)
seurat_integrated <- IntegrateLayers(object = merged_seurat, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)
seurat_integrated[["RNA"]] <- JoinLayers(seurat_integrated[["RNA"]])

# Análisis integrando los datos
seurat_integrated <- FindNeighbors(seurat_integrated, reduction = "integrated.cca", dims = 1:20)
seurat_integrated <- FindClusters(seurat_integrated, resolution = 0.5)
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:20, reduction = "integrated.cca")

DimPlot(seurat_integrated, reduction = "umap", group.by = c("Condition", "seurat_clusters"),label = T)

saveRDS(seurat_integrated,"seurat_integrated.rds")
