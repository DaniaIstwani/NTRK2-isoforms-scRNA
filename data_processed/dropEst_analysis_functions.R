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
library(rtracklayer)
library(harmony)

read_rds_to_seurat <- function(directory) {
  rds_files <- list.files(directory, pattern = "\\.rds$", full.names = TRUE) 
  
  if (length(rds_files) == 0) {
    stop("No .rds files found in the directory.")
  }
  
  print(rds_files)  # Print to check if files are detected
  
  seurat_list <- list()
  
  for (file in rds_files) {
    print(paste("Reading:", file))  # Print each file being read
    
    sample_name <- tools::file_path_sans_ext(basename(file)) 
    
    count_matrix <- readRDS(file)
    
    if (!is.list(count_matrix) || is.null(count_matrix$cm)) {
      stop(paste("File", file, "does not contain a valid count matrix."))
    }
    
    seurat_obj <- CreateSeuratObject(counts = count_matrix$cm)
    seurat_list[[sample_name]] <- seurat_obj
  }
  
  return(seurat_list)
}

rename_cells_in_seurat_objects <- function(seurat_objects_list) {
  
  #Function to extract the naming pattern (e.g., "10X05_1")
  #*input:* seurat_objects_trunc (a list of Seurat objects with names like count_matrix.rds_10X05_1.bam.1_trunc.gtf).
  #*Output:* Each Seurat object in seurat_objects_trunc now has cell names prefixed with the corresponding naming pattern (e.g., "10X05_1_Cell_1")
  
  # Extract the names of the Seurat objects
  names_list <- names(seurat_objects_list)
  
  # Extract the naming patterns using regex
  naming_patterns <- gsub("^.*_(10X\\d+_\\d+)\\..*$", "\\1", names_list)
  
  print(naming_patterns)
  
  for (i in 1:length(seurat_objects_list)) {
    prefix <- naming_patterns[i]
    
    current_cell_names <- colnames(seurat_objects_list[[i]])
    
    new_cell_names <- paste0(prefix, "_", current_cell_names)
    
    # Assign the new cell names to the Seurat object
    colnames(seurat_objects_list[[i]]) <- new_cell_names
  }
  
  # Return the modified list of Seurat objects
  return(seurat_objects_list)
}


normalize_and_merge <- function(seurat_objects_list) {
  #function that normalizes each seurat object, and merge the list into one Seurat
  #input list of seurats, output normalised and combined seurat
  # Normalize each seurat
  seurat_objects_list <- lapply(seurat_objects_list, NormalizeData)
  
  # Merge datasets
  combined_seurat <- merge(
    seurat_objects_list[[1]], 
    y = seurat_objects_list[-1], 
    merge.data = TRUE
  )
  
  return(combined_seurat)
}


update_orig_ident <- function(seurat_object) {
  #function that updates the orig.ident metadata column based on the sample names extracted from the cell names.
  #input: merged seurat object, output: updated metadata with corresponant cell names.
  # Extract sample names from cell names
  sample_names <- sub("^(10X[0-9]+_[0-9]+)_.*", "\\1", colnames(seurat_object))
  
  # Update orig.ident
  seurat_object$orig.ident <- sample_names
  
  return(seurat_object)
}

process_seurat <- function(seurat_object, 
                           nfeatures = 2000,        
                           npcs = 20,               
                           resolution = 0.5,        
                           run_umap = TRUE,         
                           umap_dims = 1:10) {      
  
  
  # Find variable features
  seurat_object <- FindVariableFeatures(seurat_object, nfeatures = nfeatures)
  
  # Scale data
  seurat_object <- ScaleData(seurat_object)
  
  # Run PCA
  seurat_object <- RunPCA(seurat_object, npcs = npcs)
  
  # Find neighbors and clusters
  seurat_object <- FindNeighbors(seurat_object, dims = 1:npcs)
  seurat_object <- FindClusters(seurat_object, resolution = resolution)
  
  # Run UMAP
  if (run_umap) {
    seurat_object <- RunUMAP(seurat_object, dims = umap_dims)
  }
  
  return(seurat_object)
}

run_harmony <- function(seurat_object) {
  # Function to integrate the Seurat object using Harmony and recompute UMAP
  # Input: Merged Seurat object
  # Output: Integrated Seurat object with UMAP computed and plotted
  
  # Run Harmony to correct for batch effects
  seurat_object <- RunHarmony(
    object = seurat_object,
    group.by.vars = "orig.ident"  # Replace "orig.ident" with the batch variable if different
  )
  
  # Recompute UMAP on the Harmony-corrected embedding
  seurat_object <- RunUMAP(
    object = seurat_object,
    reduction = "harmony",  # Use the Harmony reduction
    dims = 1:30             # Use the same number of dimensions as used in Harmony
  )
  
  # Plot the UMAP
  harmony_umap_plot <- DimPlot(
    object = seurat_object,
    reduction = "umap",
    group.by = "orig.ident"  # Color by batch (or any other metadata column)
  )
  
  print(harmony_umap_plot)  # Display the UMAP plot
  
  return(seurat_object)
}


run_CSCORE_on_Seurat_v5 <- function(seurat_object, 
                                    n_genes = 5000, 
                                    assay = "RNA", 
                                    layer = "data", 
                                    cscore_genes = NULL) {
  ##Input:
  #seurat_object: A Seurat v5 object.
  #n_genes: Number of top genes to select based on mean expression (default: 5000).
  #assay: The assay to use (default: "RNA").
  #layer: The layer to extract (default: "data" for normalized data).
  #cscore_genes: Optional. A vector of gene names to use instead of selecting top genes.
  #Output:A list containing:
  #coexpression_matrix: The CSCORE co-expression matrix.
  #p_values: The raw p-values from CSCORE.
  #adjusted_p_values: The BH-adjusted p-values.
  #selected_genes: The genes used for the analysis.
  
  # Load required libraries
  if (!require(Seurat)) {
    install.packages("Seurat")
    library(Seurat)
  }
  
  # Step 1: Extract normalized data
  normalized_data <- LayerData(seurat_object, assay = assay, layer = layer)
  normalized_data <- as.matrix(normalized_data)  # Convert to dense matrix
  
  # Step 2: Calculate mean expression
  mean_exp <- rowMeans(normalized_data)
  
  # Step 3: Select top genes
  if (is.null(cscore_genes)) {
    genes_selected <- names(sort(mean_exp, decreasing = TRUE))[1:n_genes]
  } else {
    genes_selected <- cscore_genes  # Use provided genes if available
  }
  
  # Step 4: Subset the normalized data
  normalized_data_selected <- normalized_data[genes_selected, ]
  
  # Step 5: Create a minimal Seurat object
  seurat_subset <- CreateSeuratObject(counts = normalized_data_selected, assay = assay)
  seurat_subset$nCount_RNA <- colSums(normalized_data_selected)
  
  # Step 6: Run CSCORE
  CSCORE_result <- CSCORE(seurat_subset, genes = genes_selected)
  
  # Step 7: Process CSCORE results
  CSCORE_coexp <- CSCORE_result$est
  CSCORE_p <- CSCORE_result$p_value
  
  # Adjust p-values using Benjamini-Hochberg (BH) method
  p_matrix_BH <- matrix(0, length(genes_selected), length(genes_selected))
  p_matrix_BH[upper.tri(p_matrix_BH)] <- p.adjust(CSCORE_p[upper.tri(CSCORE_p)], method = "BH")
  p_matrix_BH <- p_matrix_BH + t(p_matrix_BH)
  
  # Filter co-expression matrix based on adjusted p-values
  CSCORE_coexp[p_matrix_BH > 0.05] <- 0
  
  # Return results
  return(list(
    coexpression_matrix = CSCORE_coexp,
    p_values = CSCORE_p,
    adjusted_p_values = p_matrix_BH,
    selected_genes = genes_selected
  ))
}




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


