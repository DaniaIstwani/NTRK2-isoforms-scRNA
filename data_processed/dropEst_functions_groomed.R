## MAKE SURE TO LOAD LIBRARIES:
library(SeuratObject)
library(SeuratDisk)
library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(readr)
library(DT)
library(tidyr)
library(patchwork)
# if (!requireNamespace("rtracklayer", quietly = TRUE)) {
#     install.packages("BiocManager")
#     BiocManager::install("rtracklayer")
# }

#library(rtracklayer)



# getting a seurat object by DropEst count matrix 

count_to_seurat <- function(counts_data_path, resolution = 0.1, dims = 1:10) {
  cell_counts <- readRDS(counts_data_path) #load the count data
  seurat_obj <- CreateSeuratObject(counts = cell_counts$cm) 
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj)
  
  # Find neighbors and clusters
  seurat_obj <- FindNeighbors(seurat_obj, dims = dims)
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
  
  # Run UMAP
  seurat_obj <- RunUMAP(seurat_obj, dims = dims)
  
  # filter based on cells' mitochondrial content, feature count , RNA count
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  
  return (seurat_obj)
}

##################################################


# calculate gene experssion in a seurat
calc_gene_exp <- function(goi= gene_of_interest, seurat_obj) {
  gene_exp_df <- FetchData(seurat_obj, vars = goi)
  gene_exp_val <- sum(gene_exp_df[[goi]])
  return(gene_exp_val)
}


################################################

#plot gene expression based on seurat object:

## create a dataframe with one-gene's expression vs list of seurat objects: 


### Function to extract and aggregate gene expression data
extract_gene_expression <- function(gene_of_interest, seurat_obj_list, labels) {
  expression_data <- data.frame()  
  
  
  for (i in seq_along(seurat_obj_list)) {
    seurat_obj <- seurat_obj_list[[i]]
    if (gene_of_interest %in% rownames(seurat_obj@assays$RNA)) {
      gene_expression <- calc_gene_exp(gene_of_interest, seurat_obj)
    } else {
      # If the gene is not found, assign 0 for the dataset
      gene_expression <- 0
    }
    
    # Create a dataframe with expression values and corresponding dataset label
    temp_df <- data.frame(Expression = gene_expression, Object = labels[i])
    
    # Combine the dataframes
    expression_data <- rbind(expression_data, temp_df)
  }
  
  return(expression_data)  # Return the aggregated expression data
}

