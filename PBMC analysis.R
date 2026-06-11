# Single cell RNA

# This is the code for the analysis of the single cell RNA data from this paper:

# Machine learning approach to single cell transcriptomic analysis of Sjogren's disease reveals altered activation states of B and T lymphocytes

# This code was used for generating figures for my presentation and my report

# Date: 04/06/2026 15:15


############### Libraries ####################

library(Seurat)
library(SeuratDisk)
library(anndata)
library(dplyr)
library(readr)
library(Matrix)
library(tidyr)
library(ggplot2)

# Libraries needed for ssGSEA
library(GSVA)
library(GSEABase)
library(tidyverse)

# To add gene sets to ssGSEA results for excel sheet
library(msigdbr)


# To load in final DC object

pcDC_only <- readRDS("C:/Users/thijm/Downloads/PBMC_dc_only_data")


############### Loading in the data ###################

# Function to help loading in the data downloaded from the paper
load_patient <- function(path, donor_id) {
  
  # Load matrix
  matrix_folder_name <- paste0(donor_id, "_matrix.mtx")
  
  counts <- readMM(file.path(path, matrix_folder_name, matrix_folder_name))
  
  # Load features
  feature_folder_name <- paste0(donor_id, "_features.tsv")
  
  features <- read.delim(file.path(path, feature_folder_name, feature_folder_name),
                         header = FALSE, stringsAsFactors = FALSE)
  # Load barcodes
  barcode_folder_name <- paste0(donor_id, "_barcodes.tsv")
  file_list <- paste0(path, "/", barcode_folder_name)
  barcode_files <- list.files(file_list, full.names = FALSE)
  
  # Create new barcodes file in matrix folder to not interfere with the original barcodes
  barcode_folder_name <- paste0(donor_id, "_barcodes.tsv")
  barcode_file_path <- file.path(path, barcode_folder_name, barcode_folder_name)
  
  barcodes <- read.delim(barcode_file_path, header = FALSE)
  
  # Assign row/column names
  gene_names <- make.unique(features[,2])   
  rownames(counts) <- gene_names
  colnames(counts) <- barcodes[,1]
  
  # Create Seurat object
  raw_seurat <- CreateSeuratObject(counts = counts)
  
  # Add donor metadata
  raw_seurat$donor <- donor_id
  
  return(raw_seurat)
}


# Loading in the data

seu_PSS_1 <- load_patient("C:/Users/thijm/Downloads/GSE253568_RAW","GSM8023482_SjD_1")

seu_PSS_2 <- load_patient("C:/Users/thijm/Downloads/GSE253568_RAW","GSM8023483_SjD_2")

seu_PSS_3 <- load_patient("C:/Users/thijm/Downloads/GSE253568_RAW","GSM8023484_SjD_3")

seu_PSS_4 <- load_patient("C:/Users/thijm/Downloads/GSE253568_RAW","GSM8023488_SjD_4")

seu_PSS_5 <- load_patient("C:/Users/thijm/Downloads/GSE253568_RAW","GSM8023489_SjD_5")

seu_PSS_6 <- load_patient("C:/Users/thijm/Downloads/GSE253568_RAW","GSM8023490_SjD_6")

seu_PSS_7 <- load_patient("C:/Users/thijm/Downloads/GSE253568_RAW","GSM8023494_SjD_7")

seu_PSS_8 <- load_patient("C:/Users/thijm/Downloads/GSE253568_RAW","GSM8023495_SjD_8")

seu_PSS_9 <- load_patient("C:/Users/thijm/Downloads/GSE253568_RAW","GSM8023496_SjD_9")

# SICCA donors

seu_SICCA_1 <- load_patient("C:/Users/thijm/Downloads/GSE253568_RAW","GSM8023480_HD_1")

seu_SICCA_2 <- load_patient("C:/Users/thijm/Downloads/GSE253568_RAW","GSM8023481_HD_2")

seu_SICCA_3 <- load_patient("C:/Users/thijm/Downloads/GSE253568_RAW","GSM8023485_HD_3")

seu_SICCA_4 <- load_patient("C:/Users/thijm/Downloads/GSE253568_RAW","GSM8023486_HD_4")

seu_SICCA_5 <- load_patient("C:/Users/thijm/Downloads/GSE253568_RAW","GSM8023487_HD_5")

seu_SICCA_6 <- load_patient("C:/Users/thijm/Downloads/GSE253568_RAW","GSM8023491_HD_6")

seu_SICCA_7 <- load_patient("C:/Users/thijm/Downloads/GSE253568_RAW","GSM8023492_HD_7")

seu_SICCA_8 <- load_patient("C:/Users/thijm/Downloads/GSE253568_RAW","GSM8023493_HD_8")


# Creating a list with all donor data objects
seurat_list <- list(seu_PSS_1,seu_PSS_2,seu_PSS_3,seu_PSS_4,seu_PSS_5,seu_PSS_6,seu_PSS_7, seu_PSS_8, seu_PSS_9,
                   seu_SICCA_1, seu_SICCA_2, seu_SICCA_3, seu_SICCA_4, seu_SICCA_5, seu_SICCA_6, seu_SICCA_7, seu_SICCA_8)

# Setting donor ids
names(seurat_list) <- c("Sjogren_1", "Sjogren_2", "Sjogren_3", "Sjogren_4", "Sjogren_5", "Sjogren_6", "Sjogren_7", "Sjogren_8", "Sjogren_9",
                        "SICCA_1", "SICCA_2", "SICCA_3", "SICCA_4", "SICCA_5", "SICCA_6", "SICCA_7", "SICCA_8")

# Creating the complete seurat object
# Use seurat_list[[1]] to get a seurat object with Sjogren_1 and seurat_list[-1] to get a list with all other donors to merge into a seurat object
PBMC_paper_1_data <- merge(x = seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list))

PBMC_paper_1_data <- JoinLayers(PBMC_paper_1_data)


# Creating donor_id meta data column

# Splitting the column names
parts <- strsplit(colnames(PBMC_paper_1_data), "_")

# Grabbing the parts with the donor ids and putting it into a metadata column
PBMC_paper_1_data$donor_id <- sapply(parts, function(x) paste(x[1:2], collapse = "_"))

# adding a metadata column for condition
PBMC_paper_1_data$condition <- sapply(parts, `[`, 1)

############### Start data analysis ##################


# Checking the data quality
VlnPlot(PBMC_paper_1_data, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
# Max features = 4000 and nCount_RNA goes up to 30000 but most below 15000

FeatureScatter(PBMC_paper_1_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# QC for the data
PBMC_paper_1_data <- subset(PBMC_paper_1_data, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA < 20000)

# Normalizing the data
PBMC_paper_1_data <- NormalizeData(PBMC_paper_1_data) # normalization.method = "LogNormalize", scale.factor = 10000

# Finding features that have high cell to cell variance. Choosing 2000 features
PBMC_paper_1_data <- FindVariableFeatures(PBMC_paper_1_data, selection.method = "vst", nfeatures = 2000)

# Scaling only top 2000 variable features 
PBMC_paper_1_data <- ScaleData(PBMC_paper_1_data, features = VariableFeatures(PBMC_paper_1_data)) 
# To prevent exploding file size only scale on top features

PBMC_paper_1_data <- RunPCA(PBMC_paper_1_data, features = VariableFeatures(object = PBMC_paper_1_data))

DimPlot(PBMC_paper_1_data, reduction = "pca")
# A lot a lot of cells so Rstudio cant plot all

############### Umapping ###################

# Creating an elbowplot to see which PCs to use
ElbowPlot(PBMC_paper_1_data)
# PC 6 cutoff

# Creating a UMAP
PBMC_paper_1_data <- FindNeighbors(PBMC_paper_1_data, dims = 1:6)

PBMC_paper_1_data <- FindClusters(PBMC_paper_1_data, resolution = 0.5)

PBMC_paper_1_data <- RunUMAP(PBMC_paper_1_data, dims = 1:6)

# Plotting the UMAP
umap_plot <- DimPlot(PBMC_paper_1_data, reduction = "umap")
umap_plot

############### DC search ######################

# List with marker genes for different types of dendritic cells
marker_panel <- list(
  pDC = c("CLEC4C", "TCF4", "GZMB", "IL3RA", "CD33", "CLEC4C", "NRP1"),
  cDC1 = c("BTLA", "CADM1", "CLEC9A", "XCR1", "BATF3"),
  cDC2 = c("CLEC10A", "CD1C", "FCER1A"),
  DC3  = c("CD14", "CD163", "CD1C", "S100A8", "S100A9","FCER1A"),
  moDC = c("C5AR1", "FCAR", "CCR2", "S100A8", "S100A9", "ITGAM", "MRC1"),
  fDC = c("CR2", "CR1", "CXCL13", "FCER2", "TMEM119", "SOX9")
)

# Checking which clusters contain dendritic cells
DotPlot(PBMC_paper_1_data, features = c(unique(unlist(marker_panel)))) +
  RotatedAxis()

# Plotting the umap and featureplots to investigate which clusters contain DCs 
umap_plot + FeaturePlot(PBMC_paper_1_data, features= c("C5AR1", "FCAR", "CCR2", "S100A8", "S100A9", "ITGAM", "MRC1"))

# pDC cluster 12, tiny part away from big cluster
# cDC1 also tiny part in cluster 12
# cDC2 also clearly in cluster 12
# DC3s bit in 12 maybe also some in 6
# moDC seem like they are in 6 and 12 

# Can't find fDCs. TMEM119 and SOX9 show almost no expression, same for CXCL13. 
# The other marker genes do show some expression but they are not only specific for fDCs

# 6 is possibly maybe macrophage cluster, including it as it might contain DCs but keep in mind

# Creating a new object with only DC-containing clusters
first_subset <- subset(PBMC_paper_1_data, idents = c(6,12))

############### Analysing the first subset ##############

# 10880 cells

# Finding variable features and scaling the data
first_subset <- FindVariableFeatures(first_subset, selection.method = "vst", nfeatures = 2000)

first_subset<- ScaleData(first_subset, features = VariableFeatures(first_subset)) 

# Running PCA
first_subset <- RunPCA(first_subset, features = VariableFeatures(object = first_subset))

DimPlot(first_subset, reduction = "pca")

# Validating whether the clusters do contain DCs
DotPlot(first_subset, features = unique(unlist(marker_panel))) +
  RotatedAxis()

# Creating an elbowplot to find which PCs to use
ElbowPlot(first_subset, ndims = 40)
# PC 16 cutoff

# Creating a UMAP
first_subset <- FindNeighbors(first_subset, dims = 1:16)

first_subset <- FindClusters(first_subset, resolution = 0.5)

first_subset <- RunUMAP(first_subset, dims = 1:16)

# Plotting the UMAP
umap_plot <- DimPlot(first_subset, reduction = "umap")
umap_plot

# Saving this intermediate file for convenience
SaveSeuratRds(first_subset, file = "C:/Users/thijm/Downloads/PBMC_intermediate_pot_dc_data")

# Checking the clusters for dendritic cells
DotPlot(first_subset, features = c(unique(unlist(marker_panel)))) +
  RotatedAxis()

# Plotting the umap and featureplots to investigate which clusters contain DCs 
umap_plot + FeaturePlot(first_subset, features= c("C5AR1", "FCAR", "CCR2", "S100A8", "S100A9", "ITGAM", "MRC1","VCAN","FCER1A"))

umap_plot + FeaturePlot(first_subset, features= c("CD14", "FCGR3A", "TNFRSF1A", "TNFRSF1B", "FCGR1A","IL1RN"))

VlnPlot(first_subset, features = c("CD14", "FCGR3A", "TNFRSF1A", "TNFRSF1B", "FCGR1A","IL1RN"))
# Seems like monocytes are the issue, also in the paper you can see big monocyte cluster

# 7 does look interesting, highest peaks for FCAR and C5AR1
# pDC is cluster 13
# cDC1 in cluster 10 (but not a lot)
# cDC2 is cluster 10
# DC3 
# moDC not 10, 11 or 13. probably not 7
# fDC again don't really see them in this data set

# don't take 0 and 1 look very monocyte from the violinplots
# there are other clusters which are suspicious but will remove in next subset
# cluster 4 has lots of FCGR3A which is a monocyte marker so also exclude that

# everything except 0,1,4 and 11

second_subset <- subset(first_subset, idents = c(2,3,5,6,7,8,9,10,12,13))


############### Analysing the second subset ################

# Finding variable features and scaling the data
second_subset <- FindVariableFeatures(second_subset, selection.method = "vst", nfeatures = 2000)

second_subset <- ScaleData(second_subset, features = VariableFeatures(second_subset)) 

# Running PCA
second_subset <- RunPCA(second_subset, features = VariableFeatures(object = second_subset))

DimPlot(second_subset, reduction = "pca")

# Creating an elbowplot to find which PCs to take
ElbowPlot(second_subset, ndims = 40)
# Hard to tell which PC is the cutoff

# using computational method to determine elbowpoint

# Get variance explained per PC (as percentages)
variance_exp <- second_subset[["pca"]]@stdev / sum(second_subset[["pca"]]@stdev) * 100

# Finds the PC where the variance explained suddenly drops (looks where drop bigger than 0.1 %)
drops <- variance_exp[-length(variance_exp)] - variance_exp[-1]
# grabs last big drop (after which the plot flattens out)
elbow_pc <- sort(which(drops > 0.1), decreasing = TRUE)[1] + 1
elbow_pc

# gives PC 14 as a cutoff

# Creating a UMAP
second_subset <- FindNeighbors(second_subset, dims = 1:14)

second_subset <- FindClusters(second_subset, resolution = 0.5)

second_subset <- RunUMAP(second_subset, dims = 1:14)

# Plotting the UMAP
umap_plot <- DimPlot(second_subset, reduction = "umap")
umap_plot


# Checking which clusters contain dendritic cells
VlnPlot(second_subset, features = c("C5AR1", "FCAR", "CCR2", "S100A8", "S100A9", "ITGAM", "MRC1","CD14","FCGR1A","VCAN","FCGR3A"))

umap_plot + FeaturePlot(second_subset, features = c("CD4","CCR7","CCR6","CCR3","CCR5","CCR4")) +
  RotatedAxis()

# pDC is cluster 10 
# cDC1 in cluster 8
# cDC2 is cluster 8
# DC3 looks like a lot of clusters not 10 or 8 probably not 9 but all other clusters could contain them
# moDC same as for DC3
# fDC probably still not present


# cluster 0,2 and 5 have low FCAR and C5AR1 expression (maybe not DCs?)
# cluster 1,3,7 maybe 4 pop up with high FCAR and C5AR1 expression

# cluster 2 seems to have part NK cells ("CD3E"-,"CD3G"-,"CD3D"-,"NCAM1+")

# the other part of cluster 2 seems like CD8+ T cells so safe to remove cluster ("CD8A","CD8B","CD4")
# cluster 0 also maybe T cells? ("CCR7+")

# cluster 5 is probably B cell ("MS4A1", "CD27", "SPN", "CD70"-,"CD24")

# Lot of clusters look relatively similar
# 5 maybe non inflammatory monocytes

# not taking 0,2,5,9

# doing a findmarkers to differentiate between monocytes and moDC clusters

# Finding marker genes
second_subset.markers <- FindAllMarkers(second_subset, only.pos = TRUE)

second_subset$clusters <- Idents(second_subset)

# For finding useful markers
filtered_markers <- second_subset.markers %>%
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
  dplyr::select(cluster, rank, gene) %>%
  pivot_wider(names_from = cluster, values_from = gene)

# cluster 8 and 10 look like DCs

# Creating a new object with only DC-containing clusters
third_subset <- subset(second_subset, idents = c(1,3,4,6,7,8,10))

############### Analysing the third subset ##################

# Finding variable features and scaling the data
third_subset <- FindVariableFeatures(third_subset, selection.method = "vst", nfeatures = 2000)

third_subset <- ScaleData(third_subset, features = VariableFeatures(third_subset)) 

# Running PCA
third_subset <- RunPCA(third_subset, features = VariableFeatures(object = third_subset))

DimPlot(third_subset, reduction = "pca")

# Creating a heatmap with all cells in the object
DimHeatmap(third_subset, dims = 1, cells = 2730, balanced = TRUE)

# Creating an elbowplot used to find which PCs to use
ElbowPlot(third_subset, ndims = 40)
# Hard to determine the elbowpoint

# Using the computational method 

# Get variance explained per PC (as percentages)
variance_exp <- third_subset[["pca"]]@stdev / sum(third_subset[["pca"]]@stdev) * 100

# Finds the PC where the variance explained suddenly drops (looks where drop bigger than 0.1 %)
drops <- variance_exp[-length(variance_exp)] - variance_exp[-1]
# grabs last big drop (after which the plot flattens out)
elbow_pc <- sort(which(drops > 0.1), decreasing = TRUE)[1] + 1
elbow_pc

# elbow finder code suggest 15

# Creating a UMAP
third_subset <- FindNeighbors(third_subset, dims = 1:15)

third_subset <- FindClusters(third_subset, resolution = 0.5)

third_subset <- RunUMAP(third_subset, dims = 1:15)

# Plotting the UMAP
umap_plot <- DimPlot(third_subset, reduction = "umap")
umap_plot

# Checking which cluster contain dendritic cells
VlnPlot(third_subset, features = c("C5AR1", "FCAR", "CCR2", "S100A8", "S100A9", "ITGAM", "MRC1"))

umap_plot + FeaturePlot(third_subset, features = c("C5AR1", "FCAR", "CCR2", "S100A8", "S100A9", "ITGAM", "MRC1")) +
  RotatedAxis()

# pDC is cluster 8 
# cDC1 some cells in cluster 5
# cDC2 is cluster 5
# DC3 and moDC both again look like all clusters except 5 and 8

# cluster 1 shows low S100A8 and 9 expression also low C5AR1 and FCAR so likely not DC

# finding gene markers to aid in clusters identification
third_subset.markers <- FindAllMarkers(third_subset, only.pos = TRUE)

third_subset$clusters <- Idents(third_subset)

# For finding useful gene markers
filtered_markers <- third_subset.markers %>%
  filter(pct.1 >= 0.20, (pct.1 - pct.2) >= 0.15, avg_log2FC >= 0.25) %>%
  # Filter: expressed in ≥20% of cluster, cluster specificity >= 0.15, meaningful effect log2FC >= 0.25
  
  # Keeping genes with padj < 0.05
  filter(p_val_adj < 0.05) %>%

  # Selecting top 10 per cluster
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%        
  slice_head(n = 10) %>%
  ungroup()

# Creating a table with top 10 DE genes per cluster
top10_table <- filtered_markers %>%
  group_by(cluster) %>%
  mutate(rank = row_number()) %>%
  dplyr::select(cluster, rank, gene) %>%
  pivot_wider(names_from = cluster, values_from = gene)


# cluster 0 see S100A8 and 9 VCAN, APOBEC3A, CCR1, HLA-A, B

# cluster 2 and 4 high adj p val for C5AR1 and FCAR, likely moDC

# 9 seems like T cell CD3E, D and G very clearly

# for cluster 7 don't see any of the DC markers, cluster 6 shows lot of MT genes (don't include)
# cluster 3 do see some HLA genes but also no marker genes i would expect

# don't include cluster 1, 6, 7 and 9


############### Analysing only cDC1,2 and pDC object ##############

# As you would not expect to get crystal clear moDC clusters as they only completely differentiate in the tissue and this is PBMC data 
# I don't include those clusters in the analysis as they also contain monocytes

# Creating an object without moDC and DC3s as that annotation is not crystal clear unlike for cDC1, 2 and pDC
# Especially due to finding out that monocytes/moDC cluster per donor
pcDC_only <- subset(third_subset, idents = c(5,8))

# 328 cells in total

# Finding variable features and scaling the data
pcDC_only <- FindVariableFeatures(pcDC_only, selection.method = "vst", nfeatures = 2000)

pcDC_only <- ScaleData(pcDC_only, features = VariableFeatures(pcDC_only)) 

# Running PCA
pcDC_only <- RunPCA(pcDC_only, features = VariableFeatures(object = pcDC_only))

DimPlot(pcDC_only, reduction = "pca")

# Creating an elbowplot used to determine which PCs to use
ElbowPlot(pcDC_only, ndims = 40)
# Unclear elbowpoint so using the computational method

# Get variance explained per PC (as percentages)
variance_exp <- pbmc_dc_only[["pca"]]@stdev / sum(pbmc_dc_only[["pca"]]@stdev) * 100

# Finds where the PC where the variance explained suddenly drops (looks where drop bigger than 0.1 %)
# grabs last big drop (after which the plot flattens out)
drops <- variance_exp[-length(variance_exp)] - variance_exp[-1]
elbow_pc <- sort(which(drops > 0.1), decreasing = TRUE)[1] + 1
elbow_pc

# PC 10 given by the elbowfinder

# Creating a UMAP
pcDC_only <- FindNeighbors(pcDC_only, dims = 1:10)

pcDC_only <- FindClusters(pcDC_only, resolution = 0.5)

pcDC_only <- RunUMAP(pcDC_only, dims = 1:10)

# Plotting the UMAP
umap_plot <- DimPlot(pcDC_only, reduction = "umap")
umap_plot

# Running findmarkers to annotate the clusters
pcDC_only.markers <- FindAllMarkers(pcDC_only, only.pos = TRUE)

pcDC_only$clusters <- Idents(pcDC_only)

# For finding useful markers
filtered_markers <- pcDC_only.markers %>%
  filter(pct.1 >= 0.20, (pct.1 - pct.2) >= 0.15, avg_log2FC >= 0.25) %>%
  # Filter: expressed in ≥20% of cluster, cluster specificity >= 0.15, meaningful effect log2FC >= 0.25
  
  # Keeping genes with padj < 0.05
  filter(p_val_adj < 0.05) %>%

  # Selecting top 10 per cluster
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 10) %>%
  ungroup()

# Creating a table with top 10 DE genes per cluster
top10_table <- filtered_markers %>%
  group_by(cluster) %>%
  mutate(rank = row_number()) %>%
  dplyr::select(cluster, rank, gene) %>%
  pivot_wider(names_from = cluster, values_from = gene)

# Creating a table for the presentation
write.csv(top10_table, "PBMC_top_10.csv", row.names = TRUE)

# Creating an excel table with all significant DE genes from the findmarkers on pcDC_data 
write.csv(filtered_markers, "PBMC_filtered_markers_10.csv", row.names = TRUE)

# Checking the clusters (For all DC marker genes)
VlnPlot(pcDC_only, features = c("CLEC10A", "CD1C", "FCER1A"))

# cluster 3 is pDC
# tiny amount of cDC1 in cluster 0/1
# cluster 0 and 1 seem like cDC2, cluster 3 also expresses cDC2 markers but especially CD1C to lower extent.
# check whether cDC2 clusters are separated by condition

# Checking the clusters with the featureplot and UMAP
umap_plot + FeaturePlot(pcDC_only, features = c("BTLA", "CADM1", "CLEC9A", "XCR1", "BATF3")) +
  RotatedAxis()


# Annotating the clusters
new.cluster.ids <- c("Primed_cDC2","Homeostatic_cDC2","DC3","pDC")
names(new.cluster.ids) <- levels(pcDC_only)
pcDC_only <- RenameIdents(pcDC_only, new.cluster.ids)

# So the 2 cDC2 clusters are separated based on functional maturity/primed state. Homeostatic_cDC2 is the more homeostatic/resting state
# Showing more HLAs and many many RP genes in its marker genes. While Primed_cDC2 showed more migration and antigen presentation genes
# indicating a primed state.

# Altering the condition to change healthy to SICCA (just plain wrong by authors)

pcDC_only$condition[pcDC_only$condition == "healthy"] <- "SICCA"

# Altering the donor_ids to say SICCA_1 etc instead of healthy_1

pcDC_only$donor_id <- sub("^healthy", "SICCA", pcDC_only$donor_id)

# Saving actual final dc object (can be loaded in at the top)
SaveSeuratRds(pcDC_only, file = "C:/Users/thijm/Downloads/PBMC_dc_only_data")


############### Donor plot generation ####################

# For all 3 papers taking 0,60 and 0,015 for the axis limits for better comparison

# Getting the DC counts per donor
cell_counts <- as.list(table(pcDC_only$donor_id))

# Creating a dataframe with donor ids and DC counts
patient_cell_df <- data.frame(donor_id = names(cell_counts), count = unlist(cell_counts))

# Getting the donor_ids and the condition from the metadata
donor_info <- unique(pcDC_only@meta.data[, c("donor_id", "condition")])

# Merging the donor_ids, condition and DC_counts into one dataframe
patient_disease_df <- merge(patient_cell_df, donor_info[, c("donor_id", "condition")], by = "donor_id", all.x = TRUE)

# Plotting the DC counts per donor into a barplot, coloured by condition
patient_plot <- ggplot(patient_disease_df, aes(x = donor_id, y = count, fill = condition)) +
  labs(x = "Donor ID", y = "Dendritic cell counts") +
  ggtitle("Dendritic cell counts per donor") +
  geom_col() +
  scale_fill_manual(values = c("SICCA" = "#04BADE", "Sjogren" = "#404040")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,60)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Changing the title and axis sizes and position
plot1 <- patient_plot + theme(plot.title = element_text(size=14, face="bold.italic", hjust = 0.5),
                              axis.title.x = element_text(size=14, face="bold"),
                              axis.title.y = element_text(size=14, face="bold"))
plot1

# Determining fraction of DCs of total cell counts

# Getting total cell counts per donor (from the complete data set)
total_cell_counts <- as.list(table(PBMC_paper_1_data$donor_id))

# Creating a dataframe with donor ids and total cell counts
total_patient_cell_df <- data.frame(donor_id = names(total_cell_counts), total_count = unlist(total_cell_counts))

# Merging the donor id, condition and Dc count dataframe with the total dc count containing dataframe
DC_total_cells_df <- merge.data.frame(patient_disease_df, total_patient_cell_df)

# Creating a barplot with the total cell counts per donor, coloured by condition
total_patient_plot <- ggplot(DC_total_cells_df, aes(x = donor_id, y = total_count, fill = condition)) +
  labs(x = "Donor ID", y = "Total cell counts") +
  ggtitle("Total cell counts per donor") +
  geom_col() +
  scale_fill_manual(values = c("SICCA" = "#04BADE", "Sjogren" = "#404040")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,10000)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Changing the title and axis title sizes and position
plot2 <- total_patient_plot + theme(plot.title = element_text(size=14, face="bold.italic", hjust = 0.5),
                                    axis.title.x = element_text(size=14, face="bold"),
                                    axis.title.y = element_text(size=14, face="bold"))
plot2

# Making a barplot showing the fraction of DCs of total cells

# Calculating the fraction of DCs of total cell counts per donor
DC_total_cells_df$fraction <- DC_total_cells_df$count/DC_total_cells_df$total_count

# Creating a barplot with the fraction of DCs of total cell counts per donor, coloured by condition
DC_fraction_plot <- ggplot(DC_total_cells_df, aes(x = donor_id, y = fraction, fill = condition)) +
  labs(x = "Donor ID", y = "Fraction of DCs from total") +
  ggtitle("Fraction of DCs from total cells per donor") +
  geom_col() +
  scale_fill_manual(values = c("SICCA" = "#04BADE", "Sjogren" = "#404040")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,0.015)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Changing the title and axis title sizes and position
plot3 <- DC_fraction_plot + theme(plot.title = element_text(size=14, face="bold.italic", hjust = 0.5),
                                  axis.title.x = element_text(size=14, face="bold"),
                                  axis.title.y = element_text(size=14, face="bold"))
plot3

# To put all 3 plots together
plot1 + plot2 + plot3


############### ssGSEA ###############


# Creating the input object needed for ssGSEA

# To be sure putting the annotations into a metadata column
pcDC_only$cell_type <- Idents(pcDC_only)

# Loading the GO_BP C5 gene set
gene_sets <- getGmt("C:/Users/thijm/Downloads/c5.go.bp.v2025.1.Hs.symbols.gmt")

# Creating an expression matrix: genes x cells
mat <- as.matrix(GetAssayData(pcDC_only, slot = "data"))

# Running ssGSEA
param <- ssgseaParam(exprData = mat, geneSets = gene_sets, normalize = TRUE)

# Calculating ssGSEA scores
ssgsea_scores <- gsva(param)

# Converting ssGSEA matrix to long format
df_ssgsea <- ssgsea_scores %>%
  as.data.frame() %>%
  rownames_to_column("pathway") %>%
  pivot_longer(cols = -pathway, names_to = "cell", values_to = "score")

# Getting all metadata column from the DC-object
meta <- pcDC_only@meta.data %>%
  rownames_to_column("cell")

# Adding cell_type and condition column to dataframe
df_ssgsea <- df_ssgsea %>%
  left_join(meta[c("cell_type","cell","condition")], by = "cell")

# Run Wilcoxon test per cluster × pathway
# Exact = false to remove errors
pvals <- df_ssgsea %>%
  group_by(cell_type, pathway) %>%
  summarise(p_value = wilcox.test(score ~ condition, exact = FALSE)$p.value)

# Adding p values to the ssGSEA results dataframe
df_ssgsea <- df_ssgsea %>%
  left_join(pvals, by = c("cell_type", "pathway"))

# Getting the mean ssGSEA scores per subtype
df_cluster_condition <- df_ssgsea %>%
  group_by(cell_type, condition, pathway) %>%
  summarise(mean_score = mean(score), .groups = "drop")

# Adding p values to the condition dataframe
df_cluster_condition <- df_cluster_condition %>%
  left_join(pvals, by = c("cell_type", "pathway"))

# To show difference between Sjogren's and SICCA setting SICCA to 0
df_cluster_condition <- df_cluster_condition %>%
  group_by(cell_type, pathway) %>%
  mutate(norm_score = mean_score - mean_score[condition == "SICCA"])

# Spliting the mean and norm score into 2 columns (one for Sjogren's and one for SICCA)
df_diff <- df_cluster_condition %>%
  pivot_wider(id_cols = c(cell_type, pathway, p_value), names_from = condition, values_from = c(mean_score, norm_score))

# Keeping pathways with at least one significant result for one of the DC types
df_diff <- df_diff %>%
  group_by(pathway) %>%
  filter(any(p_value < 0.05)) %>%
  ungroup()

# Getting the top pathways 
top_pathways <- df_diff %>%
  arrange(norm_score_Sjogren) %>%
  slice_head(n = 20) %>%
  pull(pathway)
top_pathways

# To change pathways for the plot to some other interesting pathways
top_pathways <- c("GOBP_INTERFERON_ALPHA_PRODUCTION",	"GOBP_INFLAMMATORY_RESPONSE","GOBP_ADAPTIVE_IMMUNE_RESPONSE")

top_pathways <- c("GOBP_RESPONSE_TO_INTERFERON_ALPHA","GOBP_RESPONSE_TO_INTERFERON_BETA","GOBP_INTERLEUKIN_27_MEDIATED_SIGNALING_PATHWAY")

# Filters to only the top pathways
df_GO_BP_plot <- df_cluster_condition %>%
  filter(pathway %in% top_pathways)


# Creating significance labels (p val <0.05 *, p val <0.01 **, p val <0.001 *** otherwise empty string)
sig_label <- case_when(
  df_GO_BP_plot$p_value < 0.001 ~ "***",
  df_GO_BP_plot$p_value < 0.01  ~ "**",
  df_GO_BP_plot$p_value < 0.05  ~ "*",
  TRUE                          ~ ""
)

# Adding significance labels to the dataframe
df_GO_BP_plot$sig_label <- sig_label

# Only include labels for Sjogren's bars in the plot
df_GO_BP_plot$sig_label[df_GO_BP_plot$condition != "Sjogren"] <- ""

# Removes GOBP_ from the start of every pathway name
df_GO_BP_plot$pathway <- sub("^GOBP_", "", df_GO_BP_plot$pathway)

# Creating a grouped barplot showing the normalized ssGSEA scores per Sjogren's DC type
plot3 <- ggplot(df_GO_BP_plot, aes(x = cell_type, y = norm_score, fill = condition)) +
  geom_col(position = position_dodge(width = 0.9)) +
  # For adding error asterix to bars
  geom_text(aes(label = sig_label),
            position = position_dodge(width = 0.8), vjust = -0.5, size = 4) +
  # Splitting the pathways into different graphs and freeing the y axis scale.
  facet_wrap(~ pathway) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.text = element_text(size = 7)) +
  labs(x = "DC subtype", y = "Mean ssGSEA score", fill = "Condition")

plot3

# for putting the significant pathways in excel

# Getting all significant pathways
significant_pathways <- subset(df_diff, p_value < 0.05)
significant_pathways <- significant_pathways[,c(1,2,3,6,7)] # removing the mean score columns

# Adding gene sets as a list to the table

# Loading all genes
go_bp <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "BP") %>%
  dplyr::select(gs_name, gene_symbol)

# Creating a list 
go_bp_genes <- go_bp %>%
  group_by(gs_name) %>%
  summarise(genes = paste(unique(gene_symbol), collapse = ", "))

# joining the gene list on the original table (genes have to be present apparently)
# takes a while

# To prevent errors in left join subset to only genes which are in the significant pathways
go_bp_genes <- subset(go_bp_genes, gs_name %in% significant_pathways$pathway)

# Merging the genes to the pathways
pathways_gene_lists <- left_join(significant_pathways,go_bp_genes, by = c("pathway" = "gs_name"), relationship = "many-to-many")

# Sorting the pathways on p value
pathways_gene_lists_dc_paper_1 <- pathways_gene_lists %>% arrange(p_value)

# saving the object as a csv file to put in excel
write.csv(pathways_gene_lists_dc_paper_1, "ssGSEA_significant_paper_1.csv", row.names = FALSE)
