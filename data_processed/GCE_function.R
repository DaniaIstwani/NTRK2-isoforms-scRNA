### gene co-expression functions:

# function to get the co-expressed genes with gene of interest (goi)

get_co_exp <- function(seurat_obj, gene_of_interest, threshold) {
  expression_data <- as.matrix(seurat_obj@assay$RNA$counts)
  expression_data <- expression_data[rowSums(expression_data) > 0 , ] #remove genes with zero expression 
  goi_exp <- expression_data[gene_of_interest, ]
  
  correlations <- apply(expression_data, 1, function(x) cor(x, goi_exp, method = "pearson"))
  correlations <- correlations[names(correlations) != gene_of_interest] #remove self corr
  
  filtered_correlations <- correlations[abs(correlations) > threshold]
  
  co_expressed_genes <- names(filtered_correlations)
  co_expressed_correlations <- as.vector(filtered_correlations)
  
  data.frame(
    from = rep(gene_of_interest, length(co_expressed_genes)),
    to = co_expressed_genes,
    weight = co_expressed_correlations
  )
}


# plotting co-expression network 
