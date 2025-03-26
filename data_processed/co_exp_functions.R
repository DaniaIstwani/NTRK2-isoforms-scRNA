
extract_normalized_data <- function(seurat_object, assay = "RNA", layer = "data") {
  
  normalized_data <- LayerData(seurat_object, assay = assay, layer = layer)
  normalized_data <- as.matrix(normalized_data)
  return(normalized_data)
}


select_top_genes <- function(normalized_data, n_top_genes = 5000) {
  mean_exp <- rowMeans(normalized_data)
  genes_selected <- names(sort(mean_exp, decreasing = TRUE))[1:n_top_genes]
  return(genes_selected)
}

run_cscore <- function(seurat_subset, genes_selected) {
  CSCORE_result <- CSCORE(seurat_subset, genes = genes_selected)
  return(CSCORE_result)
}

process_cscore_results <- function(CSCORE_result, pval_cutoff = 0.05) {
  CSCORE_coexp <- CSCORE_result$est
  CSCORE_p <- CSCORE_result$p_value
  
  # Adjust p-values using Benjamini-Hochberg (BH) method
  p_matrix_BH <- matrix(0, length(genes_selected), length(genes_selected))
  p_matrix_BH[upper.tri(p_matrix_BH)] <- p.adjust(CSCORE_p[upper.tri(CSCORE_p)], method = "BH")
  p_matrix_BH <- p_matrix_BH + t(p_matrix_BH)
  
  # Filter co-expression matrix based on adjusted p-values
  CSCORE_coexp[p_matrix_BH > pval_cutoff] <- 0
  
  return(CSCORE_coexp)
}



correlate_with_target_gene <- function(CSCORE_coexp, target_gene, correlation_threshold = 0.5) {

  target_gene_coexp_cs <- CSCORE_coexp[target_gene, ]
  significant_coexp_genes_cs <- target_gene_coexp_cs[abs(target_gene_coexp_cs) > correlation_threshold]
  
  
  return(significant_correlations)
}



visualize_top_genes_heatmap <- function(coexp_matrix, significant_correlations, target_gene, n_top_genes = 20, my_colors = colorRampPalette(c("blue", "white", "red"))(100)) {
  # Get the top n genes
  top_genes <- significant_correlations$gene[1:n_top_genes]
  
  # Subset the co-expression matrix to include only the top genes
  top_genes_coexp <- coexp_matrix[top_genes, top_genes]
  
  # Create the heatmap
  pheatmap(top_genes_coexp, 
           color = my_colors,
           cluster_rows = TRUE, 
           cluster_cols = TRUE,
           fontsize_row = 10, 
           fontsize_col = 8,
           main = paste("Top", n_top_genes, "Genes Co-expressed with", target_gene),
           show_rownames = TRUE, 
           show_colnames = TRUE,
           cellwidth = 15, 
           cellheight = 15)
}



