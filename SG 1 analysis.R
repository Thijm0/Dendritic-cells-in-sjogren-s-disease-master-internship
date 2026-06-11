# Single cell RNA 

# This is the code for the analysis of the single cell RNA data from this paper:

# GZMK+CD8+ T cells Target A Specific Acinar Cell Type in Sjögren's Disease

# This code was used for generating figures for my presentation and my report

# Date: 04/06/2026 14:00

################## libraries  #################

library(Seurat)
library(SeuratDisk)
library(anndata)
library(dplyr)
library(readr)
library(ggpubr)
library(ggplot2)
library(Matrix)
library(tidyr)

# libraries for ssGSEA
library(GSVA)
library(GSEABase)
library(tidyverse)
library(msigdbr)


# To load in the final subsetted data
dc_only_1 <- readRDS("C:/Users/thijm/Downloads/DC_paper_1_final_data")

# To load in the correctly scaled dc_only_1 data set
dc_only_1_final <- readRDS("C:/Users/thijm/Downloads/DC_paper_1_final_scaled")

################## Loading in the data ######################

raw_data <- read_h5ad("C:/Users/thijm/Downloads/Master_internship_1_data_paper_1.h5ad")

# Creating a seurat object
paper_1_data <- CreateSeuratObject(counts = t(raw_data$X), meta.data = raw_data$obs)

# Loading in gene names and respective ensembl ids from the ensemble bioMart website
mapping <- read.table("C:/Users/thijm/Downloads/mart_export.txt", sep = "\t", header = TRUE)
 
 # Grabbing the ensembl ids from the data set
genes <- rownames(paper_1_data)

# Mapping ensembl ids to gene names
gene_symbols <- mapping$Gene.name[match(genes, mapping$Gene.stable.ID)]

# Making sure gene names don't occur twice
rownames(paper_1_data) <- make.unique(gene_symbols)

# To save the paper_1_data object for the spatial analysis
SaveSeuratRds(paper_1_data, file = "C:/Users/thijm/Downloads/paper_1_data")


################## Validating authors dendritic cell annotation. ##############################


# List of gene markers per dendritic cell type
marker_panel <- list(
  pDC = c("CLEC4C", "TCF4", "GZMB", "IL3RA", "CD33", "NRP1"),
  cDC1 = c("BTLA", "CADM1", "CLEC9A", "XCR1", "BATF3"),
  cDC2 = c("CLEC10A", "CD1C", "FCER1A"),
  DC3  = c("CD14", "CD163", "CD1C", "S100A8", "S100A9"),
  moDC = c("C5AR1", "FCAR", "CCR2", "S100A8", "S100A9", "ITGAM", "MRC1"),
  fDCs = c("CR2", "CR1", "CXCL13", "FCER2", "TMEM119", "SOX9")
)

# paper_1_data is the complete data set, first running pca and then clustering and validating with the gene marker panel

# Plotting some features of the data set used for QC
VlnPlot(paper_1_data, features = c("nFeature_RNA", "nCount_RNA", "n_genes", "total_counts","pct_counts_mt"), ncol = 5)

# Adjusting the active idents to stop it from showing up in figures
Idents(paper_1_data) <- "Cells"

# QC
paper_1_data <- subset(paper_1_data, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & total_counts < 100000)

# Normalizing the data
paper_1_data <- NormalizeData(paper_1_data, normalization.method = "LogNormalize", scale.factor = 10000)

# Finding variable features 
paper_1_data <- FindVariableFeatures(paper_1_data, selection.method = "vst", nfeatures = 2000)

# Scaleing the date
paper_1_data <- ScaleData(paper_1_data, features = VariableFeatures(paper_1_data))

# Running PCA
paper_1_data <- RunPCA(paper_1_data, features = VariableFeatures(object = paper_1_data))

# Plotting the PCA plot
DimPlot(paper_1_data, reduction = "pca")

# Creating a heatmap with all cells in the object
DimHeatmap(paper_1_data, dims = 1, cells = 71700, balanced = TRUE)

# Plotting an elbowplot to determine which PCs to use
ElbowPlot(paper_1_data)
# pc 9 cutoff

# Creating a UMAP
paper_1_data <- FindNeighbors(paper_1_data, dims = 1:9)

paper_1_data <- FindClusters(paper_1_data, resolution = 0.5)

paper_1_data <- RunUMAP(paper_1_data, dims = 1:9)

# Plotting the UMAP
DimPlot(paper_1_data, reduction = "umap")

# Checking which clusters contain DCs
DotPlot(paper_1_data, features = unique(unlist(marker_panel))) +
  RotatedAxis()

# Plotting the expression of cDC2 marker genes (did this for all DC types)
FeaturePlot(paper_1_data, features = c("CLEC10A", "CD1C", "FCER1A"), order = TRUE)
# pDCs and cDC1 seem to be in cluster 12, cDC2 super obvious cluster 12
# Maybe some fDCs in cluster 14

# Definitely cluster 12, cluster 14 contains some cDC1 so also include just to be sure
# Taking uncertain clusters as well as those will be filtered out later

# Creating a DC containing object
first_subset <- subset(paper_1_data, idents = c(12,14))

################## Analyzing the first subset ##################

# Finding variable features and scaling
first_subset <- FindVariableFeatures(first_subset, selection.method = "vst", nfeatures = 2000)
first_subset <- ScaleData(first_subset, features = VariableFeatures(first_subset)) 

# Running PCA
first_subset <- RunPCA(first_subset, features = VariableFeatures(object = first_subset))

# Plotting the PCA plot
DimPlot(first_subset, reduction = "pca")

# Validating that the clusters indeed contain dendritic cells
DotPlot(first_subset, features = unique(unlist(marker_panel))) +
  RotatedAxis()

# Plotting a heatmap of all cells in this subset
DimHeatmap(first_subset, dims = 1, cells = 2704, balanced = TRUE)

# Plotting an elbowplot to determine which PCs to use
ElbowPlot(first_subset)
# PC 9 cutoff

# Creating a UMAP
first_subset <- FindNeighbors(first_subset, dims = 1:9)

first_subset <- FindClusters(first_subset, resolution = 0.5)

first_subset <- RunUMAP(first_subset, dims = 1:9)

# Plotting the UMAP
umap_plot <- DimPlot(first_subset, reduction = "umap")
umap_plot

# Grouping by disease (just to check)
DimPlot(first_subset, reduction = "umap", group.by = "disease")

# Checking which clusters contain DCs
DotPlot(first_subset, features = unique(unlist(marker_panel))) +
  RotatedAxis()

# Plotting the complete umap and the expression of some marker genes to check which clusters contain DCs
umap_plot + FeaturePlot(first_subset, features = c("CD14", "CD163", "CD1C", "S100A8", "S100A9"), order = TRUE)

# Checked for all subtypes in which clusters they could be found
# Some pDCs in 4, also cDC1 in 4
# cDC2 many in 4 as well but also in 0 a bit and 9
# moDCs in 4 and 0
# no fDCs

# Determining gene markers for all the clusters
first_subset.markers <- FindAllMarkers(first_subset, only.pos = TRUE)
first_subset.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# Putting the clusters in a meta data column
first_subset$clusters <- Idents(first_subset)

# Finding useful markers
filtered_markers <- first_subset.markers %>%
  filter(pct.1 >= 0.20, (pct.1 - pct.2) >= 0.15, avg_log2FC >= 0.25) %>%
  # filters: expressed in ≥20% of cluster, cluster specificity and meaningful effect size
  
  # Keeping genes with padj < 0.05 
  filter(p_val_adj < 0.05) %>%
  
  # Selecting top 10 per cluster
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 10) %>%
  ungroup()

# Creating a table with the top 10 marker genes per cluster
top10_table <- filtered_markers %>%
  group_by(cluster) %>%
  mutate(rank = row_number()) %>%
  select(cluster, rank, gene) %>%
  pivot_wider(names_from = cluster, values_from = gene)

# from top10 table seems like only cluster 4 can be DCs

# Creating a DC containing object
second_subset <- subset(first_subset, idents = c(4))


################## Analyzing the second subset ###################

# Finding variable features and scaling the data
second_subset <- FindVariableFeatures(second_subset, selection.method = "vst", nfeatures = 2000)
second_subset <- ScaleData(second_subset, features = VariableFeatures(second_subset)) 

# Running PCA
second_subset <- RunPCA(second_subset, features = VariableFeatures(object = second_subset))

# Plotting the PCA plot
DimPlot(second_subset, reduction = "pca")

# Validating whether the clusters contain DCs
DotPlot(second_subset, features = unique(unlist(marker_panel))) +
  RotatedAxis()

# Creating a heatmap with all cells in this subset
DimHeatmap(second_subset, dims = 1, cells = 283, balanced = TRUE)

# Plotting the elbowplot to determine which PCs to use
ElbowPlot(second_subset)
# PC 8 cutoff

# Creating a UMAP
second_subset <- FindNeighbors(second_subset, dims = 1:8)

second_subset <- FindClusters(second_subset, resolution = 0.5)

second_subset <- RunUMAP(second_subset, dims = 1:8)

# Plotting a UMAP
umap_plot <- DimPlot(second_subset, reduction = "umap")
umap_plot

# Checking which clusters contain DCs
DotPlot(second_subset, features = unique(unlist(marker_panel))) +
  RotatedAxis()

# Plotting the complete umap and expression of some marker genes to analyse which clusters contain DCs
umap_plot + FeaturePlot(second_subset, features = c("MHC2","CD163", "MARCO","MSR1"), order = TRUE)

# cluster 3 is moDC
# cluster 0 and 1 could be cDC2
# bottom part of cluster 2 is cDC1

# Annotating the clusters and setting the active idents to these annotations
new.cluster.ids <- c("cDC2", "activated cDC2", "cDC1", "moDC")
names(new.cluster.ids) <- levels(second_subset)
second_subset <- RenameIdents(second_subset, new.cluster.ids)

# renaming the final Dc object
dc_only_1 <- second_subset

# Saving the final DC containing object use for further analysis (this is the file you can load in at the top)
SaveSeuratRds(dc_only_1, file = "C:/Users/thijm/Downloads/DC_paper_1_final_data")


################## Findmarkers on second_subset to annotate ###########################

# Finding gene markers on the DC object
second_subset.markers <- FindAllMarkers(second_subset, only.pos = TRUE)
second_subset.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# Adding a metadata column for the clusters
second_subset$clusters <- Idents(second_subset)

# Finding useful markers
filtered_markers <- second_subset.markers %>%
  filter(pct.1 >= 0.20, (pct.1 - pct.2) >= 0.15, avg_log2FC >= 0.25) %>%
  # expressed in ≥20% of cluster, cluster specificity and meaningful effect size
  
  # Keeping genes with padj < 0.05 
  filter(p_val_adj < 0.05) %>%

  # Selecting top 10 per cluster
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 10) %>%
  ungroup()

# Creating a table with the top 10 marker genes per cluster
top10_table <- filtered_markers %>%
  group_by(cluster) %>%
  mutate(rank = row_number()) %>%
  select(cluster, rank, gene) %>%
  pivot_wider(names_from = cluster, values_from = gene)


################## Analysis of the dc_only_1 object #################

# Creating a UMAP for the presentation
umap_plot <- DimPlot(dc_only_1, reduction = "umap")

# Running findmarkers to find the top marker genes for the presentation
dc_only_1.markers <- FindAllMarkers(dc_only_1, only.pos = TRUE)

# Putting the idents into a metadata column
dc_only_1$clusters <- Idents(dc_only_1)

# Finding useful markers
filtered_markers <- dc_only_1.markers %>%
  filter(pct.1 >= 0.20, (pct.1 - pct.2) >= 0.15, avg_log2FC >= 0.25) %>%
  # Expressed in ≥20% of cluster, cluster specificity and meaningful effect size
  
  # Keeping genes with padj < 0.05 
  filter(p_val_adj < 0.05) %>%
  
  # Selecting top 10 per cluster
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 10) %>%
  ungroup()

# Creating a table with the top 10 marker genes per cluster
top10_table <- as.data.frame(filtered_markers) %>%
  group_by(cluster) %>%
  mutate(rank = row_number()) %>%
  dplyr::select(cluster, rank, gene) %>% # code didn't want to use the correct select
  pivot_wider(names_from = cluster, values_from = gene)

# Saving the table to a csv file
write.csv(top10_table, "C:/Users/thijm/downloads/top10_table.csv", row.names=FALSE)

# Checking marker genes in the clusters
umap_plot + FeaturePlot(dc_only_1, features = c("CD1C", "FCER1A", "CLEC10A", "CD1E", "FCGR2B"))

# Also with violinplots
VlnPlot(dc_only_1, features = c("CD1C", "FCER1A", "CLEC10A"))

# Putting the DC types into a metadata column
dc_only_1$dc_subtype <- Idents(dc_only_1)

################## Plotting DC counts per donor ################


# Renaming the disease column to condition (makes more sense and aligns with dc_paper_2)
colnames(dc_only_1@meta.data)[colnames(dc_only_1@meta.data) == "disease"] <- "condition"

# Converting to character to replace Sjogren syndrome with Sjogren_syndrome and then converting back to factor
dc_only_1@meta.data$condition <- as.character(dc_only_1@meta.data$condition)

dc_only_1@meta.data$condition[dc_only_1@meta.data$condition == "Sjogren syndrome"] <- "Sjogren_syndrome"

dc_only_1@meta.data$condition <- factor(dc_only_1@meta.data$condition)

# Creating a metadata column with the DC type
dc_only_1$cell_type <- Idents(dc_only_1)

# Getting DC counts per donor
cell_counts <- as.list(table(dc_only_1$donor_id))

# Creating a df with disease and cell count per donor
patient_cell_df <- data.frame(donor_id = names(cell_counts), count = unlist(cell_counts))

# Getting the condition per donor
donor_info <- unique(dc_only_1@meta.data[, c("donor_id", "condition")])

# Merging the cell counts and condition 
patient_condition_df <- merge(patient_cell_df, donor_info[, c("donor_id", "condition")], by = "donor_id", all.x = TRUE)

# Reordering donor_id based on condition
patient_condition_df$donor_id <- factor(patient_condition_df$donor_id, levels = patient_condition_df$donor_id[order(patient_condition_df$condition)])

# Creating a barplot with the DC counts per donor, coloured by condition
patient_plot <- ggplot(patient_condition_df, aes(x = donor_id, y = count, fill = condition)) +
  labs(x = "Donor ID", y = "Dendritic cell counts") +
  ggtitle("Dendritic cell counts per donor") +
  geom_col() +
  scale_fill_manual(values = c("normal" = "#04BADE", "Sjogren_syndrome" = "#404040")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,60)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Changing the title and axis title sizes and position
plot1 <- patient_plot + theme(plot.title = element_text(size=14, face="bold.italic", hjust = 0.5),
                              axis.title.x = element_text(size=14, face="bold"),
                              axis.title.y = element_text(size=14, face="bold"))

plot1

# Creating total cell count per donor plot

# Creating a list with the total cell counts
total_cell_counts <- as.list(table(paper_1_data$donor_id))

# Creating a dataframe with the total cell counts per donor
total_patient_cell_df <- data.frame(donor_id = names(total_cell_counts), total_count = unlist(total_cell_counts))

# Merging the total cell counts and DC counts and donor ids + condition
DC_total_cells_df <- merge(patient_condition_df, total_patient_cell_df, by = "donor_id", all.x = TRUE)

# Reordering donor_id based on condition
DC_total_cells_df$donor_id <- factor(DC_total_cells_df$donor_id, levels = DC_total_cells_df$donor_id[order(DC_total_cells_df$condition)])

# Creating a barplot with the total cell counts per donor, coloured by condition
total_patient_plot <- ggplot(DC_total_cells_df, aes(x = donor_id, y = total_count, fill = condition)) +
  labs(x = "Donor ID", y = "Total cell counts") +
  ggtitle("Total cell counts per donor") +
  geom_col() +
  scale_fill_manual(values = c("normal" = "#04BADE", "Sjogren_syndrome" = "#404040")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,7500)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Changing the title and axis title sizes and position
plot2 <- total_patient_plot + theme(plot.title = element_text(size=14, face="bold.italic", hjust = 0.5),
                                    axis.title.x = element_text(size=14, face="bold"),
                                    axis.title.y = element_text(size=14, face="bold"))

# Making fraction DC plot

# Calculating the fraction of DCs of total cells
DC_total_cells_df$fraction <- DC_total_cells_df$count/DC_total_cells_df$total_count

# Creating a plot of the fraction of DCs of total cells, coloured by condition
DC_fraction_plot <- ggplot(DC_total_cells_df, aes(x = donor_id, y = fraction, fill = condition)) +
  labs(x = "Donor ID", y = "Fraction of DCs from total") +
  ggtitle("Fraction of DCs from total cells per donor") +
  geom_col() +
  scale_fill_manual(values = c("normal" = "#04BADE", "Sjogren_syndrome" = "#404040")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,0.015)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Changing the title and axis title sizes and position
plot3 <- DC_fraction_plot + theme(plot.title = element_text(size=14, face="bold.italic", hjust = 0.5),
                                  axis.title.x = element_text(size=14, face="bold"),
                                  axis.title.y = element_text(size=14, face="bold"))

plot3


################## Normalizing dc_only_1 to account for nCount_RNA differences ##########################


# Had weird ssGSEA results where every path was negative for Sjogren's 
# from looking at the nCount_RNA plots it seemed lower in Sjogren's then SICCA
# To remove the bias regressing out nCount_RNA and nFeature_RNA

dc_only_1_norm <- ScaleData(dc_only_1, vars.to.regress = c("nCount_RNA", "nFeature_RNA"))

# Running PCA
dc_only_1 <- FindVariableFeatures(dc_only_1, selection.method = "vst", nfeatures = 2000)

dc_only_1 <- ScaleData(dc_only_1, vars.to.regress = c("nCount_RNA", "nFeature_RNA"), features = VariableFeatures(dc_only_1))

dc_only_1 <- RunPCA(dc_only_1, features = VariableFeatures(object = dc_only_1))

DimPlot(dc_only_1, reduction = "pca")

# To remove outliers in PCA
pca <- Embeddings(dc_only_1, "pca")

# Keep cells within reasonable range for PC1/PC2
keep <- abs(pca[,1]) < 20 & abs(pca[,2]) < 20  # adjust thresholds
dc_only_1 <- subset(dc_only_1, cells = rownames(pca)[keep])

# Plotting an elbowplot to find which PCs to use
ElbowPlot(dc_only_1, ndims = 40)

# Hard to point out specific elbow so using code

# Get variance explained per PC (as percentages)
variance_exp <- dc_only_1[["pca"]]@stdev / sum(dc_only_1[["pca"]]@stdev) * 100

# Finds the PC where the variance explained suddenly drops (looks where drop bigger than 0.1 %)
# grabs last big drop (after which the plot flattens out)
drops <- variance_exp[-length(variance_exp)] - variance_exp[-1]
elbow_pc <- sort(which(drops > 0.1), decreasing = TRUE)[1] + 1
elbow_pc


# Umapping
dc_only_1 <- FindNeighbors(dc_only_1, dims = 1:15)

dc_only_1 <- FindClusters(dc_only_1, resolution = 0.5)

dc_only_1 <- RunUMAP(dc_only_1, dims = 1:15)

DimPlot(dc_only_1, reduction = "umap")

# To find DE genes per cluster
dc_only_1.markers <- FindAllMarkers(dc_only_1, only.pos = TRUE)
dc_only_1.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

dc_only_1$clusters <- Idents(dc_only_1)

# Finding top markers, but ones that are actually useful
filtered_markers <- dc_only_1.markers %>%
  filter(pct.1 >= 0.20, (pct.1 - pct.2) >= 0.15, avg_log2FC >= 0.25) %>%
  # Filters: expressed in ≥20% of cluster, cluster specificity and meaningful effect size
  
  # Keeping genes with padj < 0.05 
  filter(p_val_adj < 0.05) %>%
  
  # Selecting top 10 per cluster
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 10) %>%
  ungroup()

# Changed top10 to filtered_markers
top10_table <- filtered_markers %>%
  group_by(cluster) %>%
  mutate(rank = row_number()) %>%
  dplyr::select(cluster, rank, gene) %>%
  pivot_wider(names_from = cluster, values_from = gene)

VlnPlot(dc_only_1, features = c(unlist(unique(marker_panel))))

# Cluster 2 is moDC
# Cluster 3 is cDC1
# Cluster 0 and 1 are both cDC2

# adding annotations and setting them to the active idents
new.cluster.ids <- c("cDC2", "activated cDC2", "moDC", "cDC1")
names(new.cluster.ids) <- levels(dc_only_1)
dc_only_1 <- RenameIdents(dc_only_1, new.cluster.ids)

# The final subsetted dc object
dc_only_1_final <- dc_only_1

# To save the final object (this can be loaded in at the top)
SaveSeuratRds(dc_only_1_final, file = "C:/Users/thijm/Downloads/DC_paper_1_final_scaled")



##################### ssGSEA #########################

# Running SCT transform to account for nCount_RNA differences between Sjogren's and SICCA

dc_only_1_final <- SCTransform(dc_only_1_final, verbose = FALSE)

colnames(dc_only_1_final@meta.data)[colnames(dc_only_1_final@meta.data) == "disease"] <- "condition"

# converting to character to replace Sjogren syndrome and then converting back to factor
dc_only_1_final@meta.data$condition <- as.character(dc_only_1_final@meta.data$condition)

dc_only_1_final@meta.data$condition[dc_only_1_final@meta.data$condition == "Sjogren syndrome"] <- "Sjogren_syndrome"

# Setting the condition label correct
dc_only_1_final@meta.data$condition[dc_only_1_final@meta.data$condition == "normal"] <- "SICCA"

dc_only_1_final@meta.data$condition <- factor(dc_only_1_final@meta.data$condition)

# To get cell_type in order
dc_only_1_final$cell_type <- Idents(dc_only_1_final)

Idents(dc_only_1_final) <- dc_only_1_final$cell_type

# Loading gene sets
gene_sets <- getGmt("C:/Users/thijm/Downloads/c5.go.bp.v2025.1.Hs.symbols.gmt")

# Setting the default assay to SCT
DefaultAssay(dc_only_1_final) <- "SCT"
mat <- as.matrix(GetAssayData(dc_only_1_final, layer = "data"))

# Running ssGSEA
param <- ssgseaParam(
  exprData = mat,
  geneSets = gene_sets,
  normalize = TRUE) 

ssgsea_scores <- gsva(param)

# Converting ssGSEA matrix to long format
df_ssgsea <- ssgsea_scores %>%
  as.data.frame() %>%
  rownames_to_column("pathway") %>%
  pivot_longer(cols = -pathway, names_to = "cell", values_to = "score")

# Getting all metadata column from dc_only_1
meta <- dc_only_1@meta.data %>%
  rownames_to_column("cell")

# Adding cell_type and condition column to dataframe
df_ssgsea <- df_ssgsea %>%
  left_join(meta[c("cell_type","cell","condition")], by = "cell")

# Run Wilcoxon test per cluster × pathway
# Exact = false to remove errors
pvals <- df_ssgsea %>%
  group_by(cell_type, pathway) %>%
  summarise(p_value = wilcox.test(score ~ condition, exact = FALSE)$p.value)

# Grouping pathways by cell type and condition
df_cluster_condition <- df_ssgsea %>%
  group_by(cell_type, condition, pathway) %>%
  summarise(mean_score = mean(score), .groups = "drop")

# Adding P values
df_cluster_condition <- df_cluster_condition %>%
  left_join(pvals, by = c("cell_type", "pathway"))

# To show difference between Sjogrens and SICCA, SICCA set to 0
df_cluster_condition <- df_cluster_condition %>%
  group_by(cell_type, pathway) %>%
  mutate(norm_score = mean_score - mean_score[condition == "SICCA"])

# Splits the mean and norm score into 2 columns (one for Sjogren and one for SICCA)
df_diff <- df_cluster_condition %>%
  pivot_wider(
    id_cols = c(cell_type, pathway, p_value),
    names_from = condition,
    values_from = c(mean_score, norm_score)
  )

# Getting the top pathways 
top_pathways <- df_diff %>%
  ungroup() %>%
  arrange(desc(norm_score_Sjogren_syndrome)) %>%
  slice_head(n = 20) %>%
  pull(pathway)
top_pathways

# Grabbing significant pathways
significant_pathways_grouped <- df_diff %>%
  group_by(pathway) %>%
  arrange(p_value) %>%
  filter(p_value < 0.05)

# To manually select interesting pathways
top_pathways <- c("GOBP_INTERLEUKIN_17_MEDIATED_SIGNALING_PATHWAY", "GOBP_NEGATIVE_REGULATION_OF_RESPONSE_TO_TYPE_II_INTERFERON", "GOBP_INTERLEUKIN_27_MEDIATED_SIGNALING_PATHWAY")

top_pathways <- c("GOBP_REGULATION_OF_INTERLEUKIN_1_MEDIATED_SIGNALING_PATHWAY", "GOBP_TYPE_III_INTERFERON_PRODUCTION", "GOBP_INTERLEUKIN_1_BETA_PRODUCTION")
 
# Filtering to only top pathways
df_GO_BP_plot <- df_cluster_condition %>% 
  filter(pathway %in% top_pathways)

# Creating significance labels
sig_label <- case_when(
  df_GO_BP_plot$p_value < 0.001 ~ "***",
  df_GO_BP_plot$p_value < 0.01  ~ "**",
  df_GO_BP_plot$p_value < 0.05  ~ "*",
  TRUE                          ~ ""
)

# Adding labels
df_GO_BP_plot$sig_label <- sig_label

# Only labels for Sjogren rows
df_GO_BP_plot$sig_label[df_GO_BP_plot$condition != "Sjogren_syndrome"] <- ""

# Removes GOBP_ from the start of every pathway name
df_GO_BP_plot$pathway <- sub("^GOBP_", "", df_GO_BP_plot$pathway)

# Creating grouped barplot
plot4 <- ggplot(df_GO_BP_plot, aes(x = cell_type, y = norm_score, fill = condition)) +
  geom_col(position = position_dodge(width = 0.9)) +

  # For adding error asterix
  geom_text(aes(label = sig_label), position = position_dodge(width = 0.8), vjust = 1.5, size = 4) +

  # Splitting the pathways into different graphs and freeing the y axis scale.
  facet_wrap(~ pathway) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.text = element_text(size = 7)) +
  labs(x = "DC subtype", y = "Mean ssGSEA score", fill = "Condition")

plot4
