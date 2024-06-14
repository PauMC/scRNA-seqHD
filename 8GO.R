# Analisis Gene ontology (GO) para los tipos celulares con mas cambios
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

# "CD4"-------------------------------------------------------------------------------------------------
de_markers_linfTyNK_CD4<-read_csv("~/sc-RNAseq1/analisis/3DE/de_markers_linfTyNK_CD4_MAST.csv")
de_markers_linfTyNK_CD4<-as.data.frame(de_markers_linfTyNK_CD4)
de_markers_linfTyNK_CD4 <- de_markers_linfTyNK_CD4[, -1]
rownames(de_markers_linfTyNK_CD4)<- de_markers_linfTyNK_CD4$gene

## UPREGULATED
# Filtrar genes 
genes_filtrados <- rownames(de_markers_linfTyNK_CD4[de_markers_linfTyNK_CD4$p_val_adj <= 0.05 & de_markers_linfTyNK_CD4$avg_log2FC > 0, ])

# Leer el archivo de datos biomart
biomart <- read.csv("mart_export38.txt")

# Filtrar el biomart para los genes de interes
biomart_filtrado <- biomart[biomart$Gene.name %in% genes_filtrados, ]

# Realizar el analisis de enriquecimiento de GO
GO_results <- enrichGO(gene = unique(biomart_filtrado$Gene.stable.ID), OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "ALL")

# Escribir los resultados en un archivo
# write.table(GO_results, "linfTyNK_CD4_go.csv")

# Crear el grafico de barras para los resultados de GO
barplot(GO_results, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scale = "free") +
  ggtitle("GO_results linfTyNK_CD4") + theme(axis.text = element_text(size = 40)) +
  theme(legend.text = element_text(size = 20))
dotplot(GO_results, showCategory=10,font.size = 9, label_format = 45)+facet_wrap("ONTOLOGY")+ggtitle("GO_results linfTyNK_CD4 upregulated")


# "CD8"-------------------------------------------------------------------------------------------------
de_markers_linfTyNK_CD8<-read_csv("~/sc-RNAseq1/analisis/3DE/de_markers_linfTyNK_CD8_MAST.csv")
de_markers_linfTyNK_CD8<-as.data.frame(de_markers_linfTyNK_CD8)
de_markers_linfTyNK_CD8 <- de_markers_linfTyNK_CD8[, -1]
rownames(de_markers_linfTyNK_CD8)<- de_markers_linfTyNK_CD8$gene

## DOWNREGULATED
# Filtrar genes
genes_filtrados <- rownames(de_markers_linfTyNK_CD8[de_markers_linfTyNK_CD8$p_val_adj <= 0.05 & de_markers_linfTyNK_CD8$avg_log2FC < 0, ])

# Leer el archivo de datos biomart
biomart <- read.csv("mart_export38.txt")

# Filtrar el biomart para los genes de interes
biomart_filtrado <- biomart[biomart$Gene.name %in% genes_filtrados, ]

# Realizar el analisis de enriquecimiento de GO
GO_results <- enrichGO(gene = unique(biomart_filtrado$Gene.stable.ID), OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "ALL")

# Escribir los resultados en un archivo
# write.table(GO_results, "linfTyNK_CD8_go.csv")

# Crear el grafico de barras para los resultados de GO
barplot(GO_results, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scale = "free") +
  ggtitle("GO_results linfTyNK_CD8") + theme(axis.text = element_text(size = 40)) +
  theme(legend.text = element_text(size = 20))
dotplot(GO_results, showCategory=10,font.size = 9, label_format = 45)+facet_wrap("ONTOLOGY")+ggtitle("GO_results linfTyNK_CD8 downregulated")
