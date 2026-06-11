# Single cell RNA

# This is the code for the analysis of the single cell RNA data from this paper:

# Molecular and spatial analysis of tertiary lymphoid structures in Sjogren’s syndrome

# This code was used for generating figures for my presentation and my report

# Date: 04/06/2026 

################## libraries #################

library(Seurat)
library(SeuratDisk)
library(anndata)
library(dplyr)
library(readr)
library(Matrix)
library(tidyr)
library(ggplot2)
library(limma)

# libraries needed for ssGSEA
library(GSVA)
library(GSEABase)
library(tidyverse)
library(msigdbr)


# to read in final Dc subsetted object with PSS_B

dc_only <- readRDS("C:/Users/thijm/Downloads/DC_paper_2_PSS_B_data")


################## Loading in the data ###################

# Creating a function to load in donor files and to create seurat object for each one
load_patient <- function(path, donor_id) {
  
  # Loading matrix
  matrix_folder_name <- paste0(donor_id, "_matrix.mtx")
  
  counts <- readMM(file.path(path, matrix_folder_name, matrix_folder_name))
  
  # Loading features
  feature_folder_name <- paste0(donor_id, "_features.tsv")
  
  features <- read.delim(file.path(path, feature_folder_name, feature_folder_name),
                         header = FALSE, stringsAsFactors = FALSE)
  # Loading barcodes
  barcode_folder_name <- paste0(donor_id, "_barcodes.tsv")
  file_list <- paste0(path, "/", barcode_folder_name)
  barcode_files <- list.files(file_list, full.names = FALSE)
  
  
  # Create new barcodes file in matrix folder to not interfere with the original barcodes
  barcode_folder_name <- paste0(donor_id, "_barcodes.tsv")
  barcode_file_path <- file.path(path, barcode_folder_name, barcode_folder_name)
  
  barcodes <- read.delim(barcode_file_path, header = FALSE)
  
  # Assigning row/column names
  gene_names <- make.unique(features[,2])   
  rownames(counts) <- gene_names
  colnames(counts) <- barcodes[,1]
  
  # Creating Seurat object
  raw_seurat <- CreateSeuratObject(counts = counts)
  
  # Adding donor metadata
  raw_seurat$donor <- donor_id
  
  return(raw_seurat)
}

# Loading in the different Sjogren's donors 

seu_PSS_A <- load_patient("C:/Users/thijm/Downloads/GSE272409_RAW","GSM8401511_PSS_A")

seu_PSS_B <- load_patient("C:/Users/thijm/Downloads/GSE272409_RAW","GSM8401512_PSS_B")

seu_PSS_D <- load_patient("C:/Users/thijm/Downloads/GSE272409_RAW","GSM8401513_PSS_D")

seu_PSS_E <- load_patient("C:/Users/thijm/Downloads/GSE272409_RAW","GSM8401514_PSS_E")

seu_PSS_F <- load_patient("C:/Users/thijm/Downloads/GSE272409_RAW","GSM8401515_PSS_F")

seu_PSS_G <- load_patient("C:/Users/thijm/Downloads/GSE272409_RAW","GSM8401516_PSS_G")

seu_PSS_H <- load_patient("C:/Users/thijm/Downloads/GSE272409_RAW","GSM8401517_PSS_H")

# Loading in the different SICCA patients

seu_SICCA_A <- load_patient("C:/Users/thijm/Downloads/GSE272409_RAW","GSM8401518_SICCA_A")

seu_SICCA_B <- load_patient("C:/Users/thijm/Downloads/GSE272409_RAW","GSM8401519_SICCA_B")

seu_SICCA_D <- load_patient("C:/Users/thijm/Downloads/GSE272409_RAW","GSM8401520_SICCA_D")

seu_SICCA_F <- load_patient("C:/Users/thijm/Downloads/GSE272409_RAW","GSM8401521_SICCA_F")

seu_SICCA_G <- load_patient("C:/Users/thijm/Downloads/GSE272409_RAW","GSM8401522_SICCA_G")

seu_SICCA_H <- load_patient("C:/Users/thijm/Downloads/GSE272409_RAW","GSM8401523_SICCA_H")


# Creating a combined seurat object with all patients

seurat_list <- list(seu_PSS_A, seu_PSS_B, seu_PSS_D, seu_PSS_E, seu_PSS_F, seu_PSS_G, seu_PSS_H, seu_SICCA_A, seu_SICCA_B, seu_SICCA_D, seu_SICCA_F, seu_SICCA_G, seu_SICCA_H)

# Setting donor ids
names(seurat_list) <- c("PSS_A","PSS_B","PSS_D","PSS_E","PSS_F","PSS_G","PSS_H","SICCA_A","SICCA_B","SICCA_D","SICCA_F","SICCA_G","SICCA_H")

# Merging all seurat objects and renaming to donor ids
DC_paper_2_data <- merge(x = seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list))

DC_paper_2_data <- JoinLayers(DC_paper_2_data)

# Splitting the column names
parts <- strsplit(colnames(DC_paper_2_data), "_")

# Adding a metadata column with the donor ids
DC_paper_2_data$donor_id <- sapply(parts, function(x) paste(x[1:2], collapse = "_"))

# Adding a metadata column with the condition
DC_paper_2_data$condition <- sapply(parts, `[`, 1)

################## Data analysis #################

# looking at nFeature and nCount of the data to determine QC cutoffs
VlnPlot(DC_paper_2_data, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
    
FeatureScatter(DC_paper_2_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# QC to remove potentially dead cells or doublets
DC_paper_2_data <- subset(DC_paper_2_data, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA < 250000)

# Normalizing the data
DC_paper_2_data <- NormalizeData(DC_paper_2_data)

# Determining which features have high cell to cell variance. Choosing 2000 features
DC_paper_2_data <- FindVariableFeatures(DC_paper_2_data, selection.method = "vst", nfeatures = 2000)

# Scaling only top 2000 variable features to limit exploding object size
DC_paper_2_data <- ScaleData(DC_paper_2_data, features = VariableFeatures(DC_paper_2_data)) 

# Running PCA
DC_paper_2_data <- RunPCA(DC_paper_2_data, features = VariableFeatures(object = DC_paper_2_data))

DimPlot(DC_paper_2_data, reduction = "pca")


################## Start cluster analysis/umapping ####################

# Creating an elbowplot to determine which PCs to use
ElbowPlot(DC_paper_2_data)
# PC 9 cutoff

# Creating a UMAP
DC_paper_2_data <- FindNeighbors(DC_paper_2_data, dims = 1:9)

DC_paper_2_data <- FindClusters(DC_paper_2_data, resolution = 0.5)

DC_paper_2_data <- RunUMAP(DC_paper_2_data, dims = 1:9)

# Plotting the UMAP
umap_plot <- DimPlot(DC_paper_2_data, reduction = "umap")
umap_plot

# Annotating the clusters based on a marker panel for DCs 
marker_panel <- list(
  pDC = c("CLEC4C", "TCF4", "GZMB", "IL3RA", "CD33", "NRP1"),
  cDC1 = c("BTLA", "CADM1", "CLEC9A", "XCR1", "BATF3"),
  cDC2 = c("CLEC10A", "CD1C", "FCER1A"),
  DC3  = c("CD14", "CD163", "CD1C", "S100A8", "S100A9","FCER1A"),
  moDC = c("C5AR1", "FCAR")
)

# Checking which clusters contain DCs
DotPlot(DC_paper_2_data, features = c(unique(unlist(marker_panel)))) +
  RotatedAxis()

# Plotting the umap and featureplots to investigate which clusters contain DCs 
umap_plot + FeaturePlot(DC_paper_2_data, features= c("CLEC4C", "TCF4", "GZMB", "IL3RA", "CD33", "CLEC4C", "NRP1"))
# moDC cluster 17
# cDC1 cluster 17 and 14 (very spread out)
# cDC2 14 and 17
# pDC maybe cluster 14/17

# Subsetting the seurat object to DC containing clusters
first_subset <- subset(DC_paper_2_data, idents = c(14,17))


################## Analyzing the first potentially DC subsetted data ##############


# Finding variable features and scaling the data
first_subset <- FindVariableFeatures(first_subset, selection.method = "vst", nfeatures = 2000)

first_subset <- ScaleData(first_subset, features = VariableFeatures(first_subset)) 

# Running PCA
first_subset <- RunPCA(first_subset, features = VariableFeatures(object = first_subset))

DimPlot(first_subset, reduction = "pca")

# Validating whether the clusters really contain DCs
DotPlot(first_subset, features = unique(unlist(marker_panel))) +
  RotatedAxis()

# Plotting an elbowplot to determine which PCs to use
ElbowPlot(first_subset)
# PC 5 cutoff

# Creating a UMAP
first_subset <- FindNeighbors(first_subset, dims = 1:5)

first_subset <- FindClusters(first_subset, resolution = 0.5)

first_subset <- RunUMAP(first_subset, dims = 1:5)

# Plotting the UMAP
umap_plot <- DimPlot(first_subset, reduction = "umap")

# Checking which clusters contain DCs
DotPlot(first_subset, features = unique(unlist(marker_panel))) +
  RotatedAxis()

umap_plot + FeaturePlot(first_subset, features = c("CD14", "CD163", "CD1C", "S100A8", "S100A9"))
# cDC1 in 5 (still weird spread out)
# cDC2 5, 8 and 3
# pDC couple in cluster 8
# moDC 4 and 7 (also macrophage likely)
# DC3 cluster 4


# Creating a DC-containing object
second_subset <- subset(first_subset, idents = c(3,4,5,7,8))


################## second round of clustering ###################

# Finding variable features and scaling the data
second_subset <- FindVariableFeatures(second_subset, selection.method = "vst", nfeatures = 2000)

second_subset <- ScaleData(second_subset, features = VariableFeatures(second_subset)) 

# Running PCA
second_subset <- RunPCA(second_subset, features = VariableFeatures(object = second_subset))

DimPlot(second_subset, reduction = "pca")

# Creating an elbowplot to determine which PCs to use
ElbowPlot(second_subset)
# PC 5 cutoff

# Creatign a UMAP
second_subset <- FindNeighbors(second_subset, dims = 1:5)

second_subset <- FindClusters(second_subset, resolution = 0.5)

second_subset <- RunUMAP(second_subset, dims = 1:5)

# Plotting the UMAP
umap_plot <- DimPlot(second_subset, reduction = "umap")

# Renaming the object
DC_data <- second_subset


################## Another round of clustering #############

# Checking which clusters contain DCs
DotPlot(DC_data, features = unique(unlist(marker_panel))) +
  RotatedAxis()

umap_plot + FeaturePlot(DC_data, features = c("CD14", "CD163", "CD1C", "S100A8", "S100A9"))
# cDC2 2 and part of 0
# cDC1 part of 5
# moDC 1 part of 2 and 3


# Hard to determine so will use DE genes
DC_data.markers <- FindAllMarkers(DC_data, only.pos = TRUE)
DC_data.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

DC_data$clusters <- Idents(DC_data)

# For finding useful markers

filtered_markers <- DC_data.markers %>%
  filter(pct.1 >= 0.20, (pct.1 - pct.2) >= 0.15, avg_log2FC >= 0.25) %>%
  # filters: expressed in ≥20% of cluster, cluster specificity and meaningful effect size
  
  # Keeping genes with padj < 0.05 
  filter(p_val_adj < 0.05) %>%
  
  # Selecting top 10 per cluster
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%        
  slice_head(n = 10) %>%
  ungroup()

# Creating a table with top 10 marker genes
top10_table <- filtered_markers %>%
  group_by(cluster) %>%
  mutate(rank = row_number()) %>%
  dplyr::select(cluster, rank, gene) %>%
  pivot_wider(names_from = cluster, values_from = gene)

# cluster 1 looks DC3 like
# cluster 2 looks like cDC2
# cluster 3 maybe moDC
# cluster 5 has cDC1

# Wrting to csv for creating an excel table
write.csv(top10_table, "C:/Users/thijm/downloads/top10_table_temp.csv", row.names=FALSE)

# Creating a DC_containing object
Final_dc_data <- subset(DC_data, idents = c(1,2,3,5))


################## Last time clustering #################

# Finding variable features and scaling the data
Final_dc_data <- FindVariableFeatures(Final_dc_data, selection.method = "vst", nfeatures = 2000)

Final_dc_data <- ScaleData(Final_dc_data, features = VariableFeatures(Final_dc_data)) 

# Running PCA
Final_dc_data <- RunPCA(Final_dc_data, features = VariableFeatures(object = Final_dc_data))

DimPlot(Final_dc_data, reduction = "pca")

# To remove some outliers seen in the PCA plot
pca <- Embeddings(Final_dc_data, "pca")

keep <- abs(pca[,1]) < 25 & abs(pca[,2]) < 25

Final_dc_data <- subset(Final_dc_data, cells = rownames(pca)[keep])

# Creating a heatmap with all cell from the object
DimHeatmap(Final_dc_data, dims = 1, cells = 450, balanced = TRUE)

# Creating an elbowplot to determine which PCs to sue
ElbowPlot(Final_dc_data)
# PC 13 cutoff

# Creating a UMAP
Final_dc_data <- FindNeighbors(Final_dc_data, dims = 1:13)

Final_dc_data <- FindClusters(Final_dc_data, resolution = 0.5)

Final_dc_data <- RunUMAP(Final_dc_data, dims = 1:13)

# Plotting the UMAP
umap_plot <- DimPlot(Final_dc_data, reduction = "umap")

################## Running findmarkers for easier identification ################

# Finding gene markers
Final_dc_data.markers <- FindAllMarkers(Final_dc_data, only.pos = TRUE)
Final_dc_data.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)


Final_dc_data$clusters <- Idents(Final_dc_data)

# For finding useful markers
filtered_markers <- Final_dc_data.markers %>%
  filter(pct.1 >= 0.20, (pct.1 - pct.2) >= 0.15, avg_log2FC >= 0.25) %>%
  # filters: expressed in ≥20% of cluster, cluster specificity and meaningful effect size
  
  # Keeping genes with padj < 0.05
  filter(p_val_adj < 0.05 ) %>%
  
  # Selecting top 10 per cluster
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%        
  slice_head(n = 10) %>%
  ungroup()

# Creating a table with top 10 marker genes
top10_table <- filtered_markers %>%
  group_by(cluster) %>%
  mutate(rank = row_number()) %>%
  dplyr::select(cluster, rank, gene) %>%
  pivot_wider(names_from = cluster, values_from = gene)

# Writing to a csv for excel
write.csv(top10_table, "C:/Users/thijm/downloads/top10_table_temp_2.csv", row.names=FALSE)

# Checking which clusters contain DCs
umap_plot + FeaturePlot(Final_dc_data, features = c(unique(unlist(marker_panel))))


# 1 is activated cDC2, 2 is cDC2, 5 is cDC1. (0 is macrophages, 3 has a lot of MT genes and 4 is T cells)
dc_B_only <- subset(Final_dc_data, idents = c(1,2,5))


################## dc_B_only reclustering ########################


# Finding variable features and scaling the data
dc_B_only <- FindVariableFeatures(dc_B_only, selection.method = "vst", nfeatures = 2000)

dc_B_only <- ScaleData(dc_B_only, features = VariableFeatures(dc_B_only)) 

# Running PCA
dc_B_only <- RunPCA(dc_B_only, features = VariableFeatures(object = dc_B_only))

DimPlot(dc_B_only, reduction = "pca")

# Creating an elbowplot to determine which PCs to use
ElbowPlot(dc_B_only)
# PC 6 cutoff

# Creating a UMAP
dc_B_only <- FindNeighbors(dc_B_only, dims = 1:6)

dc_B_only <- FindClusters(dc_B_only, resolution = 0.5)

dc_B_only <- RunUMAP(dc_B_only, dims = 1:6)

# Plotting the UMAP
umap_plot <- DimPlot(dc_B_only, reduction = "umap")
umap_plot

# Checking which clusters contain DCs
VlnPlot(dc_B_only, features = c(unique(unlist(marker_panel)),"APOE", "C1QA"))

# cluster 4 is cDC1
# cluster 0 cDC2
# 3 looks the most like DC3 or moDC
# cluster 1 and 2 look like macrophage clusters
# no real pDC cluster i think

# Validating the cluster annotation by finding marker genes
dc_B_only.markers <- FindAllMarkers(dc_B_only, only.pos = TRUE)
dc_B_only.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# For finding usefull gene markers
filtered_markers <- dc_B_only.markers %>%
  filter(pct.1 >= 0.20, (pct.1 - pct.2) >= 0.15, avg_log2FC >= 0.25) %>%
  # filters: expressed in ≥20% of cluster, cluster specificity and meaningful effect size
  
  # Keeping genes with padj < 0.05 
  filter(p_val_adj < 0.05) %>%
  
  # Selecting top 10 per cluster
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%        
  slice_head(n = 10) %>%
  ungroup()

# Creating a table with the top 10 marker genes
top10_table <- filtered_markers %>%
  group_by(cluster) %>%
  mutate(rank = row_number()) %>%
  dplyr::select(cluster, rank, gene) %>%
  pivot_wider(names_from = cluster, values_from = gene)


# 0 is cDC2 
# 4 is cDC1
# 3 is moDC

# 1 and 2 look like macrophages

# Creating a DC containing object
dc_B_only <- subset(dc_B_only, idents = c(0,3,4))
 
# Annotating the clusters
new.cluster.ids <- c("cDC2","moDC", "cDC1")
names(new.cluster.ids) <- levels(dc_B_only)
dc_B_only <- RenameIdents(dc_B_only, new.cluster.ids)

# Adding a metadta column with the DC types
dc_B_only$dc_subtype <- Idents(dc_B_only)

# Saving the final DC containing object. This is the object used for further analysis
SaveSeuratRds(dc_B_only, file = "C:/Users/thijm/Downloads/DC_paper_2_PSS_B_data")

# Writing to csv to create a table for the presentation
write.csv(top10_table, "C:/Users/thijm/downloads/top10_table.csv", row.names=FALSE)


################## Actual DC analysis ######################

# The dc_only has 3 clusters with cDC1, cDC2 and mono-DC.

# renaming the dc_B_only object to dc_only
dc_only <- dc_B_only


# creating a UMAP for the presentation
umap_plot <- DimPlot(dc_only, reduction = "umap")
umap_plot

# Doing findmarkers to find the top marker genes for all clusters also for the presentation
dc_only.markers <- FindAllMarkers(dc_only, only.pos = TRUE)
dc_only.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# Adding a metadata column with the clusters of dc_only
dc_only$clusters <- Idents(dc_only)

# Finding top markers
filtered_markers <- dc_only.markers %>%
  filter(pct.1 >= 0.20, (pct.1 - pct.2) >= 0.15, avg_log2FC >= 0.25) %>%
  # expressed in ≥20% of cluster, cluster specificity and meaningful effect size
  
  # Keeping genes with padj < 0.05 
  filter(p_val_adj < 0.05) %>%
  
  # Selecting top 10 per cluster
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 10) %>%
  ungroup()

# Changed top10 to filtered_markers
top10_table <- as.data.frame(filtered_markers) %>%
  group_by(cluster) %>%
  mutate(rank = row_number()) %>%
  dplyr::select(cluster, rank, gene) %>% 
  pivot_wider(names_from = cluster, values_from = gene)


# Writing to a csv for an excel table
write.csv(top10_table, "C:/Users/thijm/downloads/top10_table_dc_2.csv", row.names=FALSE)

################## Donor analysis #################


# Plotting DC counts per donor

# Creating a list with the cell counts per donor
cell_counts <- as.list(table(dc_only$donor_id))

# Creating a dataframe with the cell counts and donor ids
patient_cell_df <- data.frame(donor_id = names(cell_counts), count = unlist(cell_counts))

# Creating a dataframe with the donor ids and condition
donor_info <- unique(dc_only@meta.data[, c("donor_id", "condition")])

# Merging the donor ids condition and cell counts
patient_disease_df <- merge(patient_cell_df, donor_info[, c("donor_id", "condition")],
                            by = "donor_id", all.x = TRUE)

# Creating a plot with the cell counts per donor coloured by condition
patient_plot <- ggplot(patient_disease_df, aes(x = donor_id, y = count, fill = condition)) +
  labs(x = "Donor ID", y = "Dendritic cell counts") +
  ggtitle("Dendritic cell counts per donor") +
  geom_col() +
  scale_fill_manual(values = c("SICCA" = "#04BADE", "PSS" = "#404040")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,60)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Changing the title and axis title sizes and position
plot1 <- patient_plot + theme(plot.title = element_text(size=14, face="bold.italic", hjust = 0.5),
                     axis.title.x = element_text(size=14, face="bold"),
                     axis.title.y = element_text(size=14, face="bold"))

plot1

# Determining fraction of DCs of total cell counts

# Getting total cell counts
total_cell_counts <- as.list(table(DC_paper_2_data$donor_id))

# Creating a dataframe witht the total cell counts and donor ids
total_patient_cell_df <- data.frame(donor_id = names(total_cell_counts), total_count = unlist(total_cell_counts))

# Merging the total cell counts df with the DC counts df
DC_total_cells_df <- merge.data.frame(patient_disease_df, total_patient_cell_df)

# Creating a plot witht the total cell counts per donor, coloured by condition
total_patient_plot <- ggplot(DC_total_cells_df, aes(x = donor_id, y = total_count, fill = condition)) +
  labs(x = "Donor ID", y = "Total cell counts") +
  ggtitle("Total cell counts per donor") +
  geom_col() +
  scale_fill_manual(values = c("SICCA" = "#04BADE", "PSS" = "#404040")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,7500)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Changing the title and axis title sizes and position
plot2 <- total_patient_plot + theme(plot.title = element_text(size=14, face="bold.italic", hjust = 0.5),
                           axis.title.x = element_text(size=14, face="bold"),
                           axis.title.y = element_text(size=14, face="bold"))

# Making fraction DC plot

# Calculating the fraction of DCs from the total cells per donor
DC_total_cells_df$fraction <- DC_total_cells_df$count/DC_total_cells_df$total_count

# Creating a plot with the fraction of DCs of total cells per donor, coloured by condition
DC_fraction_plot <- ggplot(DC_total_cells_df, aes(x = donor_id, y = fraction, fill = condition)) +
  labs(x = "Donor ID", y = "Fraction of DCs from total") +
  ggtitle("Fraction of DCs from total cells per donor") +
  geom_col() +
  scale_fill_manual(values = c("SICCA" = "#04BADE", "PSS" = "#404040")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,0.015)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Changing the title and axis title sizes and position
plot3 <- DC_fraction_plot + theme(plot.title = element_text(size=14, face="bold.italic", hjust = 0.5),
                         axis.title.x = element_text(size=14, face="bold"),
                         axis.title.y = element_text(size=14, face="bold"))

plot3

################## ssGSEA2 #######################

# Creating the input object needed for ssGSEA

# Putting the Dc types into a metadata column
dc_only$cell_type <- Idents(dc_only)

# Loading the GO BP C5 gene set
gene_sets <- getGmt("C:/Users/thijm/Downloads/c5.go.bp.v2025.1.Hs.symbols.gmt")

# Creating an expression matrix: genes x cells
mat <- as.matrix(GetAssayData(dc_only, layer = "data"))

# Running ssGSEA
param <- ssgseaParam(exprData = mat, geneSets = gene_sets, normalize = TRUE)

# Calculating the ssGSEA scores
ssgsea_scores <- gsva(param)

# Converting ssGSEA matrix to long format
df_ssgsea <- ssgsea_scores %>%
  as.data.frame() %>%
  rownames_to_column("pathway") %>%
  pivot_longer(cols = -pathway, names_to = "cell", values_to = "score")

# Getting all metadata column from dc_only
meta <- dc_only@meta.data %>%
  rownames_to_column("cell")

# Adding cell_type and condition column to dataframe
df_ssgsea <- df_ssgsea %>%
  left_join(meta[c("cell_type","cell","condition")], by = "cell")

# Running Wilcoxon test per cluster × pathway
# Exact = false to remove errors
pvals <- df_ssgsea %>%
  group_by(cell_type, pathway) %>%
  summarise(p_value = wilcox.test(score ~ condition, exact = FALSE)$p.value)

# Adding p values to the dataframe
df_ssgsea <- df_ssgsea %>%
  left_join(pvals, by = c("cell_type", "pathway"))

# Getting the mean ssgsea scores per subtype
df_cluster_condition <- df_ssgsea %>%
  group_by(cell_type, condition, pathway) %>%
  summarise(mean_score = mean(score), .groups = "drop")

# Adding p values to the dataframe
df_cluster_condition <- df_cluster_condition %>%
  left_join(pvals, by = c("cell_type", "pathway"))

# To show difference between PSS and SICCA setting SICCA to 0
df_cluster_condition <- df_cluster_condition %>%
  group_by(cell_type, pathway) %>%
  mutate(norm_score = mean_score - mean_score[condition == "SICCA"])

# Splitting the mean and norm scores into 2 columns (one for PSS and one for SICCA)
df_diff <- df_cluster_condition %>%
  pivot_wider(
    id_cols = c(cell_type, pathway, p_value),
    names_from = condition,
    values_from = c(mean_score, norm_score)
  )

# Keeping pathways with at least one significant result for one of the DC types
df_diff <- df_diff %>%
  group_by(pathway) %>%
  filter(any(p_value < 0.05)) %>%
  ungroup()

# Getting the top pathways 
top_pathways <- df_diff %>%
  # ungroup() %>%
  arrange(norm_score_PSS) %>%
  slice_head(n = 20) %>%
  pull(pathway)
top_pathways

# To change the top pathways for plotting potentially interesting pathways
top_pathways <- c("GOBP_INTERFERON_ALPHA_PRODUCTION",	"GOBP_INFLAMMATORY_RESPONSE","GOBP_ADAPTIVE_IMMUNE_RESPONSE")

top_pathways <- c("GOBP_INTERLEUKIN_17_MEDIATED_SIGNALING_PATHWAY","GOBP_NEGATIVE_REGULATION_OF_RESPONSE_TO_TYPE_II_INTERFERON",
"GOBP_INTERLEUKIN_27_MEDIATED_SIGNALING_PATHWAY")


# Filtering to only the top pathways
df_GO_BP_plot <- df_cluster_condition %>%
  filter(pathway %in% top_pathways)


# Creating significance labels
sig_label <- case_when(
  df_GO_BP_plot$p_value < 0.001 ~ "***",
  df_GO_BP_plot$p_value < 0.01  ~ "**",
  df_GO_BP_plot$p_value < 0.05  ~ "*",
  TRUE                          ~ ""
)

# Adding significance labels
df_GO_BP_plot$sig_label <- sig_label

# Only labeling PSS rows
df_GO_BP_plot$sig_label[df_GO_BP_plot$condition != "PSS"] <- ""

# Removes GOBP_ from the start of every pathway name
df_GO_BP_plot$pathway <- sub("^GOBP_", "", df_GO_BP_plot$pathway)


# Creating a grouped barplot showing the normalized ssGSEA scores per celltype for Sjogren's DCs
plot3 <- ggplot(df_GO_BP_plot, aes(x = cell_type, y = norm_score, fill = condition)) +
  geom_col(position = position_dodge(width = 0.9)) +
  # For adding error asterix

  geom_text(aes(label = sig_label),
            position = position_dodge(width = 0.8), vjust = -0.5, size = 4) +
  # Splitting the pathways into different graphs and freeing the y axis scale.

  facet_wrap(~ pathway) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.text = element_text(size = 7)) +
  labs(x = "DC subtype", y = "Mean ssGSEA score", fill = "Condition")

plot3
