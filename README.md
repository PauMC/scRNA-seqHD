# scRNA-seqHD
In the quest for biomarkers for Huntington’s disease (HD), high-throughput and unbiased screening of RNA molecules can help identify cost-effective and sensitive indicators. Although the expression of mutant huntingtin occurs in blood cells, previous transcriptomics studies using peripheral blood from HD mutation carriers have yielded inconsistent candidates. This inconsistency may be due to the altered gene expression patterns being obscured by the highly heterogeneous and dynamic nature of blood. Single-cell gene profiling might offer insights into the most relevant blood cell types for biomarkers.

We conducted a single-cell RNA-seq assay on peripheral blood mononuclear cells from two groups: symptomatic HD patients and control donors matched for gender and age, following Chromium 10x Genomics procedures. To perform the bioinformatics analysis, the usual mapping workflow is followed using the Cell Ranger count software developed by 10x Genomics and the human reference genome GRCh38 optimized for scRNA-seq. The mapping was done separately for the sequencing results of the control and hd conditions. The following steps are carried out in the RStudio environment using the Seurat package, designed for quality control, analysis, and data exploration. The script used in RStudio can be found in this repository, divided into:

- calidadQC.R: Quality analysis and preprocessing
- clustering_all.R: Differentiation of clusters into monocytes and dendritic cells; T lymphocytes and NK cells; B lymphocytes
- clustering_monoDC.R: Characterization of monocyte and dendritic cell clusters
- clustering_linfTandNK.R: Characterization of T lymphocyte and NK cell clusters
- clustering_linfB.R: Characterization of B lymphocyte clusters
- DE.R: Differential expression analysis
- visualizresultsDE.R: Visualization and interpretation of the differential expression analysis results

