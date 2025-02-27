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

library(rtracklayer)

# Seurat preprocessing: 

preprocess_seurat <- function(seurat_obj, dims = 1:10, resolution = 0.1) {
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj)
  seurat_obj <- FindNeighbors(seurat_obj, dims = dims)
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
  seurat_obj <- RunUMAP(seurat_obj, dims = dims)
  
  # Calculate mitochondrial percentage
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  # Filter cells based on mitochondrial content, feature count, and RNA count
  seurat_obj <- subset(
    seurat_obj,
    subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5
  )
  
  # Return the processed Seurat object
  return(seurat_obj)
}


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

 # calculate gene experssion in a seurat
 calc_gene_exp <- function(goi= goi, seurat_obj) {
   gene_exp_df <- FetchData(seurat_obj, vars = goi)
   gene_exp_val <- sum(gene_exp_df[[goi]])
   return(gene_exp_val)
}


# visualize gene expression between different seurat objects.

# Function to extract and aggregate gene expression data
extract_gene_expression <- function(gene_of_interest, seurat_obj_list, labels) {
  expression_data <- data.frame()  # Initialize empty dataframe
  
  # Loop through each Seurat object and fetch the gene expression
  for (i in seq_along(seurat_obj_list)) {
    seurat_obj <- seurat_obj_list[[i]]
    if (gene_of_interest %in% rownames(seurat_obj@assays$RNA)) {
      gene_data <- FetchData(seurat_obj, vars = gene_of_interest)
      gene_expression <- sum(gene_data[[gene_of_interest]])  # get raw expression values for each cell
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


# Function to extract and sum gene expression for multiple genes
extract_gene_expression_multi_obj <- function(genes_of_interest, seurat_obj_list, labels) {
  expression_data <- data.frame()  
  
  for (i in seq_along(seurat_obj_list)) {
    seurat_obj <- seurat_obj_list[[i]]  
    total_expression <- 0  
    
    # Loop through each gene of interest and sum their expression
    for (gene in genes_of_interest) {
      if (gene %in% rownames(seurat_obj@assays$RNA)) {
        gene_data <- FetchData(seurat_obj, vars = gene)
        total_expression <- total_expression + sum(gene_data[[gene]])  # Add to total expression
      }
    }
    
    # Create a dataframe with total expression values and corresponding object label
      temp_df <- data.frame(Expression = total_expression, Object = labels[i], Gene = gene)
    
    # Combine the dataframes
      expression_data <- rbind(expression_data, temp_df)
  }
  
  return(expression_data)  
}
  

extract_gene_expression_multi_genes <- function(genes_of_interest, seurat_obj_list, labels) {
  # Initialize an empty list to store temporary data frames
  expression_list <- list()
  
  for (i in seq_along(seurat_obj_list)) {
    seurat_obj <- seurat_obj_list[[i]]
    
    for (gene in genes_of_interest) {
      # Reset total_expression for each gene
      total_expression <- 0 
      
      if (gene %in% rownames(seurat_obj@assays$RNA)) {
        total_expression <- sum(FetchData(seurat_obj, vars = gene), na.rm = TRUE)
        
        # Add a data frame for each gene's expression to the list
        expression_list[[length(expression_list) + 1]] <- data.frame(
          Expression = total_expression,
          Object = labels[i],
          Gene = gene  # Include gene name here
        )
      } else {
        message(paste("Warning:", gene, "not found in Seurat object", labels[i]))
      }
    }
  }
  
  # Combine all data frames into one
  expression_data <- do.call(rbind, expression_list)
  
  return(expression_data)
}


# Function to plot aggregated gene expression data
plot_gene_expression <- function(expression_data, gene_of_interest, sample_name, labels) {
  p <- ggplot(expression_data, aes(x = Object, y = Expression, fill = Object)) +
    geom_bar(stat = "identity", width = 0.3) + 
    geom_text(aes(label = round(Expression, 2)), vjust = -0.5, size = 5) + 
    theme_minimal() +
    labs(title = paste("Total Expression of", gene_of_interest, "Across", sample_name, "Object(s)"), 
         y = "Total Expression Level", x = "Object") +
    scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.1))) +  
    theme(axis.text.x = element_text(angle = 35, hjust = 1),  # Rotate x-axis labels for better readability
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.3),
          plot.margin = margin(25, 10, 10, 10)) +
    scale_fill_brewer(palette = "Set2")  #color palette for bars

  plot_file_name <- paste0(gene_of_interest, sample_name, "_barplot.png")
  #print(plot_file_name)
  
  #ggsave(plot_file_name, plot = p, width = 8, height = 7, dpi = 300, device = "png")
  
  return(p)  
}

# Function to plot aggregated gene expression data for multiple genes
# plot_gene_expression_multi <- function(expression_data, sample_name) {
#   if (!"Gene" %in% colnames(expression_data)) {
#     stop("Error: 'Gene' column not found in expression_data.")
#   }
#   
#   # Plot with ggplot
#   p <- ggplot(expression_data, aes(x = Object, y = Expression, fill = Gene)) +
#     geom_bar(stat = "identity", position = position_dodge(), width = 0.6) + 
#     geom_text(aes(label = round(Expression, 2)), 
#               position = position_dodge(0.6), 
#               vjust = -0.5, size = 4) +
#     theme_minimal() +
#     labs(title = paste("Total Expression Across", sample_name, "Object(s)"), 
#          y = "Total Expression Level", x = "Object") +
#     scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.1))) +
#     theme(axis.text.x = element_text(angle = 35, hjust = 1),
#           axis.text = element_text(size = 12),
#           axis.title = element_text(size = 14),
#           plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
#           plot.margin = margin(30, 10, 10, 10)) +
#     scale_fill_brewer(palette = "Set2")  # color palette for bars
#   
#   plot_file_name <- paste0(sample_name, "_multi_gene_expression_barplot.png")
#   ggsave(plot_file_name, plot = p)
#   
#   return(p)  
# }


# Function to rename clusters, subset, and generate violin plots
process_seurat_and_plot <- function(seurat_obj, new_cluster_ids, genes_to_plot, plot_name_prefix, sample_name) {
  
  # 1: Rename clusters based on new cluster IDs
  names(new_cluster_ids) <- levels(seurat_obj)
  seurat_obj <- RenameIdents(seurat_obj, new_cluster_ids)
  
  # 2: Subset Seurat object for the specified cell types
  subset_seurat <- subset(seurat_obj, idents = c("Astrocytes", "Neurons", "Oligo"))
  
  # 3: Generate the violin plot for the specified genes
  p <- VlnPlot(subset_seurat, features = genes_to_plot) +
    ggtitle(paste(genes_to_plot, "expression level in sample", sample_name, "-", sample_name))+
    theme(plot.title = element_text(size = 12, hjust = 0.3, face = "bold"),  
          plot.margin = margin(30, 10, 10, 10),  
          plot.title.position = "plot")
  
  # 4: Save the plot with a custom name
  plot_file_name <- paste0(plot_name_prefix, "_VlnPlot.png")
  ggsave(plot_file_name, plot = p)
  print(p)
  
  # 5: Return the Seurat object and the plot for further use
  return(list(seurat_obj = seurat_obj, plot = p))
}


# find markers for a given Seurat object
find_markers_for_object <- function(seurat_obj) {
  markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  return(markers)
}




# Function to calculate co-expression with a gene for a given cell type
get_coexpression_with_gene <- function(seurat_obj, gene_of_interest, cell_type) {
  # Subset cells of the specified cell type
  subset_obj <- subset(seurat_obj, idents = cell_type)
  
  # normalized expression matrix and convert to dense matrix
  expr_matrix <- as.matrix(GetAssayData(subset_obj, slot = "data"))
  
  # gene of interest exists?
  if (!(gene_of_interest %in% rownames(expr_matrix))) {
    stop(paste("Gene", gene_of_interest, "not found in the dataset."))
  }
  
  # get the gene vector
  gene_vector <- expr_matrix[gene_of_interest, ]
  
  # correlations
  coexpression <- cor(t(expr_matrix), gene_vector, method = "pearson") # can i get the r value?
  
  # Convert to a tidy data frame
  coexpression_df <- data.frame(
    Gene = rownames(expr_matrix),
    Correlation = coexpression
  )
  
  # Sort by correlation
  coexpression_df <- coexpression_df[order(-coexpression_df$Correlation), ]
  
  return(coexpression_df)
}


