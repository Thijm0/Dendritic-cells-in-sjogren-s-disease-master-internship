# Spatial transcriptomics

# This is the code used for the analysis of the spatial transcriptomics data of this paper:

# GZMK+CD8+ T cells Target A Specific Acinar Cell Type in Sjögren's Disease

# This code was used for generating figures for my presentation and my report

# Date: 04/06/2026 13:16

########### Libraries ##################

library(readr)
library(dplyr)
library(ggplot2)
library(Seurat)
library(anndataR) 
library(SeuratDisk)
library(SeuratObject)
library(SPOTlight)
library(SingleCellExperiment)


# To load in Visium starting seurat object (after generating it in this code)
Spatial_SG_1 <- readRDS("C:/Users/thijm/Downloads/Spatial_SG_1.rds")


########### loading in the Visium data #################

# Loading in table with gene names and ensembl codes for correcting gene annotation
mapping <- read.table("C:/Users/thijm/Downloads/mart_export(1).txt", sep = "\t", header = TRUE)

# Function which creates seurat objects from the h5ad files and corrects the counts and data layers
load_data <- function(path,file_name) {
  
  # Creating a complete path
  path_folder_name <- paste0(path, file_name)
  
  # Creating a seurat object
  seurat_object <- read_h5ad(path_folder_name, as = "Seurat")
  
  # Grabbing the counts
  raw <- seurat_object@assays$RNA@layers$X
  
  # Grabbing the barcodes and gene names from the meta dat
  barcodes <- rownames(seurat_object@meta.data)
  genes <- rownames(seurat_object)
  
  # Ensemble to gene names 
  gene_symbols <- mapping$Gene.name[match(genes, mapping$Gene.stable.ID)]
  
  # To deal with NAs replacing them with their ensemble codes again
  gene_symbols[is.na(gene_symbols)] <- genes[is.na(gene_symbols)]
  
  # Assigning the barcodes and gene names to the counts
  colnames(raw) <- barcodes
  rownames(raw) <- make.unique(gene_symbols)
  
  # Recreating the RNA assay so it matches what seurat expects
  seurat_object[["RNA"]] <- CreateAssayObject(counts = raw)
  
  return(seurat_object)
}

block_1 <- load_data("C:/Users/thijm/Downloads/","56c59f1f-aff3-46c3-8f7a-a8659973492a.h5ad")
block_2 <- load_data("C:/Users/thijm/Downloads/","e0cc7d54-59eb-456c-8210-fa285747604e.h5ad")
block_3 <- load_data("C:/Users/thijm/Downloads/","691e41cc-f866-4fdf-afde-6be35424d0ad.h5ad")
block_4 <- load_data("C:/Users/thijm/Downloads/","57273732-8149-400d-ab2a-41b7876c6ffe.h5ad")
block_5 <- load_data("C:/Users/thijm/Downloads/","bd6655ab-1acc-48bb-ad66-a72d24aecbbe.h5ad")
block_6 <- load_data("C:/Users/thijm/Downloads/","ea9a14aa-7b33-4b11-83db-273264d135de.h5ad")
block_7 <- load_data("C:/Users/thijm/Downloads/","b5e5b9f2-1d74-4590-8dd1-6fc505ca6f5b.h5ad")
block_8 <- load_data("C:/Users/thijm/Downloads/","480b1cb3-bafd-4ecb-b191-9eb06e68fd90.h5ad")
block_9 <- load_data("C:/Users/thijm/Downloads/","8cfe9833-7d20-430d-880c-a75ae14a6f40.h5ad")
block_10 <- load_data("C:/Users/thijm/Downloads/","62e1b767-90e6-4b36-9a16-6951faaba4a3.h5ad")
block_13 <- load_data("C:/Users/thijm/Downloads/","abd579ae-12e4-4629-a9c8-e248df00b41d.h5ad")
block_14 <- load_data("C:/Users/thijm/Downloads/","45ab7b40-1476-4fd8-9de3-e94c0020f1cd.h5ad")
block_15 <- load_data("C:/Users/thijm/Downloads/","7c8a030e-598b-4979-b1b2-c74a2522f205.h5ad")
block_16 <- load_data("C:/Users/thijm/Downloads/","1190d582-0403-46f9-9e41-c13a622b7add.h5ad")
block_17 <- load_data("C:/Users/thijm/Downloads/","e2e2cf7a-c99e-421e-b633-1d4561be2b44.h5ad")
block_18 <- load_data("C:/Users/thijm/Downloads/","a96bbfb0-b996-4a88-9dfa-b683be8ebc9e.h5ad")
block_21 <- load_data("C:/Users/thijm/Downloads/","a3866e01-15d0-439a-9f9a-b72640d6dadd.h5ad")
block_22 <- load_data("C:/Users/thijm/Downloads/","e65f325e-fd3e-42fc-86dc-92770b7ef4d3.h5ad")
block_24 <- load_data("C:/Users/thijm/Downloads/","f0ca2a42-8fcf-4575-ab95-c501db82ec13.h5ad")
block_25 <- load_data("C:/Users/thijm/Downloads/","7e1bf179-9486-49d9-a898-1ec20b716b7f.h5ad")
block_26 <- load_data("C:/Users/thijm/Downloads/","f3f2ab0a-94a0-4014-b0c0-834c74cadc65.h5ad")
block_27 <- load_data("C:/Users/thijm/Downloads/","198cc97b-2550-4a30-b708-cdcfc9cb6c90.h5ad")
block_30 <- load_data("C:/Users/thijm/Downloads/","85f5f467-c06c-40e5-9cb3-c189499cb967.h5ad")
block_31 <- load_data("C:/Users/thijm/Downloads/","204c83ea-76c4-4727-a8b0-1625748927f3.h5ad")


# Creating a list of all blocks for merging (needed later)
block_list <- list(block_1, block_2, block_3, block_4, block_5, block_6, block_7, block_8, block_9, block_10, block_13, block_14, block_15, block_16, block_17, block_18,
  block_21, block_22, block_24, block_25, block_26, block_27, block_30, block_31)

# Adding names to the list to be able to separate in the final object
names(block_list) <- c("block_1", "block_2", "block_3", "block_4", "block_5", "block_6", "block_7", "block_8", "block_9", "block_10", "block_13", "block_14",
                        "block_15", "block_16", "block_17", "block_18", "block_21", "block_22", "block_24", "block_25", "block_26", "block_27", "block_30", "block_31")

# Merging all seurat objects into one plus merging the spatial reduction as well
Spatial_SG_1 <- merge(x = block_list[[1]], y = block_list[-1], merge.dr = TRUE, add.cell.ids = names(block_list))

########### QC ####################

# Checking nCount and nFeatures
VlnPlot(Spatial_SG_1, features = c("nCount_RNA", "nFeature_RNA"), pt.size = 0.1)

# Decided not to subset the data to not remove any real cells. Seems like the authors have already done based on looking at their annotation.

# Normalizing the data
Spatial_SG_1 <- NormalizeData(Spatial_SG_1)

# Finding variable features
Spatial_SG_1 <- FindVariableFeatures(Spatial_SG_1, selection.method = "vst", nfeatures = 2000)

# Scaling the data
Spatial_SG_1 <- ScaleData(Spatial_SG_1) 

# Running PCA
Spatial_SG_1 <- RunPCA(Spatial_SG_1) 

# Plotting PCA
DimPlot(Spatial_SG_1, reduction = "pca")

# Umapping

ElbowPlot(Spatial_SG_1, ndims = 50)

# PC 9 cutoff 

Spatial_SG_1 <- FindNeighbors(Spatial_SG_1, dims = 1:9)

Spatial_SG_1 <- FindClusters(Spatial_SG_1)

Spatial_SG_1 <- RunUMAP(Spatial_SG_1, dims = 1:9)

# Plotting the UMAP
DimPlot(Spatial_SG_1, reduction = "umap")

# Putting the spatial coords into the meta data
coords <- Spatial_SG_1@reductions$spatial@cell.embeddings

Spatial_SG_1$imagecol <- coords[, "spatial_1"]
Spatial_SG_1$imagerow <- coords[, "spatial_2"]

# To save the final seurat object to load in
SaveSeuratRds(Spatial_SG_1, file = "C:/Users/thijm/Downloads/Spatial_SG_1.rds")
# This is what is loaded in at the start

########### Using spotlight for spot deconvolution ##############

# Creating an object which contains the cell annotations from dc_only_1_final in the complete dc_paper_1 object
# Need to load in paper_1_data and dc_only_1_final from the dc_paper_1 analysis


# To load in dc_only_1_final
dc_only_1_final <- readRDS("C:/Users/thijm/Downloads/DC_paper_1_final_scaled")

# To load in paper_1_data 
paper_1_data <- readRDS("C:/Users/thijm/Downloads/paper_1_data")

# Putting the dc subtypes in a meta data column and grabbing it
dc_only_1_final$dc_subtype <- Idents(dc_only_1_final)
dc_meta <- dc_only_1_final@meta.data["dc_subtype"]

# Adding the meta data column to the complete paper_1_data object
paper_1_data <- AddMetaData(object = paper_1_data, metadata = dc_meta)

# Converting both meta data columns to characters
paper_1_data$dc_subtype <- as.character(paper_1_data$dc_subtype)
paper_1_data$cell_type <- as.character(paper_1_data$cell_type)

# Combing the annotations with priority to the DC annotations. If the annotation is na use the authors annotation
paper_1_data$annotation_combined <- ifelse(!is.na(paper_1_data$dc_subtype), paper_1_data$dc_subtype, paper_1_data$cell_type)

# Filling counts and data of paper_1_data for findmarkers
paper_1_data <- NormalizeData(paper_1_data)
paper_1_data <- FindVariableFeatures(paper_1_data)
paper_1_data <- ScaleData(paper_1_data)

# Setting active idents to the correct annotations
Idents(paper_1_data) <- paper_1_data$annotation_combined

# Running findmarkers for spotlight. Takes quite a while without presto.
markers <- FindAllMarkers(paper_1_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)

# To save the object
saveRDS(object = markers, file = "C:/Users/thijm/Downloads/Markers_spotlight.rds")

# To load it in
markers <- readRDS("C:/Users/thijm/Downloads/Markers_spotlight.rds")
  
# Getting the top 50 markers per cell type
markers_top <- markers %>%
  group_by(cluster) %>%
  top_n(50, wt = avg_log2FC)

# Altering both paper_1_data and the spatial data to single cell experiment object for input into spotlight

paper_1_SCE <- as.SingleCellExperiment(paper_1_data)

# Settings column labels (doesn't do this automatically)
colLabels(paper_1_SCE) <- paper_1_data$annotation_combined

Spatial_SG_SCE <- as.SingleCellExperiment(Spatial_SG_1)

# Running spotlight
res <- SPOTlight(
  x = paper_1_SCE,
  y = Spatial_SG_SCE,
  mgs = markers_top,
  weight_id = "avg_log2FC",
  group_id = "cluster",
  gene_id = "gene")


# Getting matrix with proportion estimates
props <- res$mat

# Getting empty spots (without cells labelled as unknown)
empty_spots <- rownames(Spatial_SG_1@meta.data)[Spatial_SG_1$cell_type == "unknown"]

# Making spots without cells 0 (To not get estimates on spots without cells)
props[empty_spots, ] <- NA

# Setting a new metadata column for the estimated cell counts per spot.
Spatial_SG_1$estimated_cells <- NA

# Removes NA so that no cell counts are estimated for spots without cells
valid_blocks <- setdiff(unique(Spatial_SG_1$Block), NA)

# For loop that calculates per block the estimated cell count per spot based on an average of 7 cells per spot.
# Uses the nCount_RNA median per slide to determine which spots have more or less cells than this average.
for (slide in valid_blocks) {
  # Grabs all spots from a specific block
  idx <- which(Spatial_SG_1$Block == slide)
  
  # Only get spots which have high enough nCount_RNA (remove background spots)
  tissue_idx <- idx[Spatial_SG_1$nCount_RNA[idx] > 1000]
  
  # Calculate median nCount_RNA per block
  med  <- median(Spatial_SG_1$nCount_RNA[tissue_idx])
  
  # Calculate number of cells per spot per block using 7 as average
  Spatial_SG_1$estimated_cells[tissue_idx] <- (Spatial_SG_1$nCount_RNA[tissue_idx] / med) * 7
}

# Normalizing proportions to stop cell counts exceeding cap
props_norm <- props / rowSums(props)

# Multiplying cell proportions with estimated cell counts to get estimate of cell counts per spot for all cell types
cell_counts <- round(props_norm * Spatial_SG_1$estimated_cells)

# Total cell counts per spot
Spatial_SG_1$cell_counts <- rowSums(cell_counts)

# Setting a hard cap of 15 cells per spot
Spatial_SG_1$cell_counts <- pmin(rowSums(cell_counts), 15)

# Adding the cell counts per cell type as metadata columns
Spatial_SG_1 <- AddMetaData(Spatial_SG_1, cell_counts)


########### Image generation ###############

# Creating separate objects per block from the Spatial_SG_1 object

# Creating an integer for the loop
x <- 1 

# Grabbing the block IDs from the list names
block_ids <- as.numeric(gsub("block_", "", names(block_list)))

# For loop which creates subsets for each block from the Spatial_SG_1 object
for (i in seq_along(block_list)) {
  # Setting end to grab correct slice
  end <- x + 4991
  
  # Subsetting Spatial_SG_1 to the correct cells
  final <- Spatial_SG_1[,x:end]
  
  # Changing Nas to 0 for plotting purposes (spots without cells)
  final$Block[is.na(final$Block)] <- "0"
  
  # Assigning the object the correct name
  new_name <- paste0("block_",block_ids[i], "_final")
  assign(new_name, final)
  
  # Adding 4992 for the next slice/object
  x = x + 4992
}

# Next some examples of how i generated the plots, 
# of course if you want to look at different cell types you only need to change the colour to the cell type of you choice
# if you want to look at different slides just change the block name

# Used geom_point(size = 1.7) if there was a lot of overlap between spots, better visability for images with a lot of cells
# Used a bigger size so the dots slightly overlap which is easier on the eyes.

# Looking at cell proportion estimates per cell type (Setting props for points with no cells to NA allows for grey background)
ggplot(block_31_final@meta.data, aes(x = imagecol, y = imagerow, colour = cell_counts)) +
  geom_point(size = 2, alpha = 0.6) +
  coord_fixed() +
  xlim(0,6000) +
  ylim(0,6000) +
  scale_colour_viridis_c(option = "plasma") +
  xlab("X") +
  ylab("Y")
  
# Looking at cell counts for specific cell types per spot
ggplot(block_31_final@meta.data, aes(x = imagecol, y = imagerow, colour = cDC2)) +
  geom_point(size = 2, alpha = 0.6) +
  coord_fixed() +
  xlim(0,6000) +
  ylim(0,6000) +
  scale_color_continuous(palette = "plasma") +
  xlab("X") +
  ylab("Y")

# If you want to plot specific T cell types to proper names for plotting purposes
colnames(block_2_final@meta.data)[colnames(block_2_final@meta.data) == "CD4-positive,_alpha-beta_T_cell"] <- "CD4_T_cell"

# Plotting CD4_T_cells (adding breaks stops legend showing half values)
ggplot(block_2_final@meta.data, aes(x = imagecol, y = imagerow, colour = CD4_T_cell)) +
  geom_point(size = 1.7, alpha = 0.6) +
  coord_fixed() +
  xlim(0,3500) +
  ylim(0,3500) +
  scale_colour_viridis_c(option = "plasma", breaks = c(0,1,2,3,4,5,6,7,8,9,10)) +
  xlab("X") +
  ylab("Y")
