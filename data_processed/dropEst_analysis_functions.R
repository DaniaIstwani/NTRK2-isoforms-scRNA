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
    
    seurat_obj <- CreateSeuratObject(
      counts = count_matrix$cm,
      min.cells = 3,
      min.features = 200
      )
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
  # Function to normalize each Seurat object with LogNormalize and merge them
  # Input: list of Seurat objects
  # Output: merged Seurat object with LogNormalized RNA assay
  
  # Normalize each Seurat object
  seurat_objects_list <- lapply(seurat_objects_list, function(x) {
    NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
  })
  
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
  #input: merged seurat object, output: updated metadata with correspondent cell names.
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

run_harmony <- function(seurat_object, group.by.vars = "orig.ident", npcs = 30) {
  # Function to integrate the Seurat object using Harmony and recompute UMAP
  # Input: Merged Seurat object
  # Output: Integrated Seurat object with UMAP computed and plotted
  #Reference: Korsunsky et al., 2019 (Nature Methods) (why Harmony)
  
  # Run Harmony to correct for batch effects
  seurat_object <- RunHarmony(
    object = seurat_object,
    group.by.vars = group.by.vars,  
    reduction.use = "pca",
    dims = 1:npcs,
    assay.use = "RNA",  
    project.dim = FALSE,
    verbose = TRUE
  )
  
  # Use harmony reduction for UMAP
  seurat_object <- RunUMAP(
    object = seurat_object,
    reduction = "harmony",
    dims = 1:npcs,
    verbose = TRUE
  )
  
  return(seurat_object)
}


run_CSCORE_on_Seurat_v5 <- function(
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
    seurat_object, 
    n_genes = 3000,          
    assay = "RNA", 
    layer = "data", 
    cscore_genes = NULL
) {
  # Validate inputs
  if (!assay %in% names(seurat_object@assays)) {
    stop("Assay not found in Seurat object")
  }
  
  # 1. Extract sparse data (no conversion to dense)
  sparse_data <- LayerData(seurat_object, assay = assay, layer = layer)
  
  # 2. Select genes (sparse-compatible)
  if (is.null(cscore_genes)) {
    # Calculate mean expression without dense conversion
    gene_means <- Matrix::rowMeans(sparse_data)
    genes_selected <- names(sort(gene_means, decreasing = TRUE)[1:n_genes])
  } else {
    genes_selected <- intersect(cscore_genes, rownames(sparse_data))
  }
  
  # 3. Subset sparse matrix (memory-efficient)
  sparse_subset <- sparse_data[genes_selected, ]
  
  # 4. Convert ONLY the subset to dense (minimal memory impact)
  dense_subset <- as.matrix(sparse_subset)
  
  # 5. Create minimal Seurat object
  seurat_subset <- CreateSeuratObject(
    counts = dense_subset,
    assay = "RNA"
  )
  
  # 6. Run CS-CORE
  CSCORE_result <- CSCORE::CSCORE(
    object = seurat_subset,
    genes = genes_selected
  )
  
  # 7. Process results (unchanged)
  CSCORE_coexp <- CSCORE_result$est
  CSCORE_p <- CSCORE_result$p_value
  
  p_matrix_BH <- matrix(0, length(genes_selected), length(genes_selected))
  p_matrix_BH[upper.tri(p_matrix_BH)] <- p.adjust(CSCORE_p[upper.tri(CSCORE_p)], method = "BH")
  p_matrix_BH <- p_matrix_BH + t(p_matrix_BH)
  
  CSCORE_coexp[p_matrix_BH > 0.05] <- 0
  
  return(list(
    coexpression_matrix = CSCORE_coexp,
    p_values = CSCORE_p,
    adjusted_p_values = p_matrix_BH,
    selected_genes = genes_selected
  ))
}

get_high_correlation_features <- function(matrix, variable_vector, threshold, method = 'pearson'){
  # taking a matrix and a vector (numerical), find features (rows) that have absolute correlation score (both positive and negative) greater than your set threshold value with your variable vector. 
  correlations <- lapply(1:nrow(matrix), function(x) {cor(as.vector(as.numeric(matrix[x,])),as.vector(variable_vector), method=method)})
  
  correlations <- as.matrix(unlist(correlations))
  
  rownames(correlations) <- rownames(as.matrix(matrix))
  
  high_corr_features <- correlations[which(abs(correlations) > threshold), drop = FALSE]
  
  return(high_corr_features)
  
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


