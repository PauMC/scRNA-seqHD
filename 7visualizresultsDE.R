# Librerias
library(readr)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(gridExtra)

# Lista de archivos
archivos <- c(
  "~/sc-RNAseq1/analisis/3DE/de_markers_seuratB_naive_MAST.csv",
  "~/sc-RNAseq1/analisis/3DE/de_markers_seuratB_memory_MAST.csv",
  "~/sc-RNAseq1/analisis/3DE/de_markers_seuratB_intermediate_MAST.csv",
  "~/sc-RNAseq1/analisis/3DE/de_markers_monoDC_CD14_MAST.csv",
  "~/sc-RNAseq1/analisis/3DE/de_markers_monoDC_CD16_MAST.csv",
  "~/sc-RNAseq1/analisis/3DE/de_markers_monoDC_cDC2_MAST.csv",
  "~/sc-RNAseq1/analisis/3DE/de_markers_linfTyNK_CD4_MAST.csv",
  "~/sc-RNAseq1/analisis/3DE/de_markers_linfTyNK_CD8_MAST.csv",
  "~/sc-RNAseq1/analisis/3DE/de_markers_linfTyNK_CD8_naive_MAST.csv",
  "~/sc-RNAseq1/analisis/3DE/de_markers_linfTyNK_MAIT_MAST.csv",
  "~/sc-RNAseq1/analisis/3DE/de_markers_linfTyNK_NK_MAST.csv"
)

# Función para filtrar y contar filas
filtrar_y_contar <- function(archivo) {
  datos <- read_csv(archivo)
  filtrado_down <- filter(datos, p_val_adj <= 0.05 & avg_log2FC < 0)
  filtrado_up <- filter(datos, p_val_adj <= 0.05 & avg_log2FC > 0)
  return(data.frame(Tipocelular = basename(archivo), Tipo = c("Down", "Up"), Valor = c(nrow(filtrado_down), nrow(filtrado_up))))
}

# Aplicar la función a cada archivo y almacenar resultados
resultados <- lapply(archivos, filtrar_y_contar)

# Combinar resultados en un data frame
resultados_df <- do.call(rbind, resultados)
resultados_df$Tipocelular <- substr(resultados_df$Tipocelular, 12, nchar(resultados_df$Tipocelular) - 9)

# Convertir el tipo a factor para el orden correcto en el gráfico
resultados_df$Tipo <- factor(resultados_df$Tipo, levels = c("Down", "Up"))

# Gráfico de barras agrupadas
ggplot(resultados_df, aes(x = Valor, y = Tipocelular, fill = Tipo)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Número de genes upregulados y downregulados significativamente",
       x = "Número de genes",
       y = "Tipo celular") +
  scale_fill_manual(values = c("blue1", "red1"), labels = c("Down", "Up")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


# EnhancedVolcano y Gene ontology (GO) para los tipos celulares con mas cambios
# Función para crear gráficos EnhancedVolcano
createEnhancedVolcano <- function(data, title) {
  EnhancedVolcano(data, 
                  lab = data$gene,
                  x = 'avg_log2FC',
                  y = 'p_val_adj',
                  title = title,
                  xlim = c(min(data$avg_log2FC, na.rm = TRUE) - 0.5, max(data$avg_log2FC, na.rm = TRUE) + 0.5),
                  ylim = c(0, max(-log10(data$p_val_adj), na.rm = TRUE) + 0.5),
                  xlab = bquote(~Log[2] ~ "fold change"),
                  pCutoff = 0.05,
                  FCcutoff = 1,
                  cutoffLineType = 'twodash',
                  cutoffLineWidth = 0.8,
                  pointSize = 2.0,
                  labSize = 3,
                  colAlpha = 1,
                  legendLabels = c('Not sig.','Log (base 2) FC','p-adj', 'p-adj & Log (base 2) FC'),
                  legendPosition = 'right',
                  legendLabSize = 10,
                  legendIconSize = 4.0)
}

#Gene ontology (GO) 
# Biomart
biomart <- read.csv("mart_export38.txt")

# Función para realizar el análisis de enriquecimiento de GO y generar los gráficos
perform_GO_analysis <- function(data, cell_type) {
  # Filtrar los genes DE
  data_filtrado_down <- filter(data, p_val_adj <= 0.05 & avg_log2FC < 0)
  data_filtrado_up <- filter(data, p_val_adj <= 0.05 & avg_log2FC > 0)
  
  # Filtrar el biomart para los genes de interés
  biomart_filtrado_down <- biomart[biomart$Gene.name %in% data_filtrado_down$gene, ]
  biomart_filtrado_up <- biomart[biomart$Gene.name %in% data_filtrado_up$gene, ]
  
  # Realizar el análisis de enriquecimiento de GO
  GO_upregulated <- enrichGO(gene = unique(biomart_filtrado_down$Gene.stable.ID), OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "ALL")
  GO_downregulated <- enrichGO(gene = unique(biomart_filtrado_up$Gene.stable.ID), OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "ALL")
  
  # Configurar colores
  options(enrichplot.colours = c("orange", "green"))
  
  # Crear los gráficos
  a <- dotplot(GO_upregulated, showCategory = 10, font.size = 9, label_format = 45) +
    facet_wrap("ONTOLOGY") +
    ggtitle(paste("GO_results", cell_type, "upregulated"))
  
  b <- dotplot(GO_downregulated, showCategory = 10, font.size = 9, label_format = 45) +
    facet_wrap("ONTOLOGY") +
    ggtitle(paste("GO_results", cell_type, "downregulated"))
  
  return(grid.arrange(a, b, ncol = 2))
}
  
linfTyNK_CD8<-read_csv("~/sc-RNAseq1/analisis/3DE/de_markers_linfTyNK_CD8_MAST.csv")
linfTyNK_CD4<-read_csv("~/sc-RNAseq1/analisis/3DE/de_markers_linfTyNK_CD4_MAST.csv")

createEnhancedVolcano(linfTyNK_CD8, "linfTyNK_CD8")
createEnhancedVolcano(linfTyNK_CD4, "linfTyNK_CD4")

perform_GO_analysis(linfTyNK_CD8, "linfTyNK_CD8")
perform_GO_analysis(linfTyNK_CD4, "linfTyNK_CD4")
