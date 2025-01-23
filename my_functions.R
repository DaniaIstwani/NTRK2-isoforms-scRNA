

# Filtering based on minimum cell count

## Function to calculate overall proportions of each cell type in the entire dataset

calculate_overall_proportions <- function(data) {
  filtered_data <- data@meta.data %>%
    filter(Class %in% c("Astrocytes", "Neurons", "Oligos"))
  
  total_cell_count <- filtered_data %>%
    summarise(TotalCellCount = n()) %>%
    pull(TotalCellCount)
  
  cell_type_counts <- filtered_data %>%   # Calculate cell counts for each cell type
    group_by(Class) %>%
    summarise(CellCount = n(), .groups = 'drop')
  
  cell_type_counts_prop <- cell_type_counts %>%   # Calculate proportions for each cell type
    mutate(Proportion = CellCount / total_cell_count)
  
  return(cell_type_counts_prop)
}

## Function to get the proportion for a specific cell type
get_proportion <- function(proportions_df, cell_type) {
  proportion <- proportions_df %>%
    filter(Class == cell_type) %>%
    pull(Proportion)
  
  return(proportion)
}


## Function to filter cells by minimum count according to overall proportions
filter_samples_by_min_count <- function(data, min_total_count) {
  
  overall_proportions <- calculate_overall_proportions(data)
  
  # Scale minimum counts for each cell type based on overall proportions
  min_astrocyte_count <- min_total_count * get_proportion(overall_proportions, "Astrocytes")
  min_neuron_count <- min_total_count * get_proportion(overall_proportions, "Neurons")
  min_oligo_count <- min_total_count * get_proportion(overall_proportions, "Oligos")
  
  # Initialize a list to store SampleIDs that meet the criteria for each cell type
  valid_sample_ids <- list()
  
  # Check each cell type
  for (cell_type in cell_types_of_interest) {
    min_count <- switch(cell_type,
                        "Astrocytes" = min_astrocyte_count,
                        "Neurons" = min_neuron_count,
                        "Oligos" = min_oligo_count)
    
    ids <- data@meta.data %>%
      filter(Class == cell_type) %>%
      group_by(SampleID) %>%
      summarise(CellCount = n(), .groups = 'drop') %>%
      filter(CellCount >= min_count) %>%
      pull(SampleID)
    
    if (length(ids) == 0) {
      cat(sprintf("No samples have the defined minimum count for %s.\n", cell_type))
    } else {
      valid_sample_ids[[cell_type]] <- ids
    }
  }
  
  # Check if there are any SampleIDs that meet the criteria for all specified cell types
  if (length(valid_sample_ids) < 3) {
    cat("Some cell types did not meet the minimum cell count criteria.\n")
    return(NULL)
  }
  
  # Find common SampleIDs across all cell types
  samples_to_keep <- Reduce(intersect, valid_sample_ids)
  
  # Check if any samples meet the criteria
  if (length(samples_to_keep) == 0) {
    cat("No samples meet the minimum cell count across all specified cell types.\n")
    return(NULL)
  }
  
  # Subset the Seurat object to include only samples that meet the criteria
  seurat_subset <- subset(data, subset = SampleID %in% samples_to_keep & Class %in% cell_types_of_interest)
  
  # Verify if subset is empty
  if (seurat_subset@meta.data %>% nrow() == 0) {
    cat("The resulting subset of the data has no cells. Check the filtering criteria.\n")
    return(NULL)
  }
  
  # Cell Type per sample
  target_cell_counts_by_sample <- table(seurat_subset@meta.data$SampleID, seurat_subset@meta.data$Class)
  print(target_cell_counts_by_sample)
  
  # Convert the table to a data frame for plotting
  target_counts_df <- as.data.frame(target_cell_counts_by_sample)
  names(target_counts_df) <- c("Sample", "CellType", "Count")
  target_counts_df <- target_counts_df %>%
    filter(CellType %in% c("Astrocytes", "Oligos", "Neurons"))
  
  return(samples_to_keep)
}


## function to plot the kept samples vs the cell count of each type

plot_cell_type_counts <- function(counts_df) {
  ggplot(counts_df, aes(x = Sample, y = Count, fill = CellType)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Cell Type Counts by Sample",
         x = "Sample",
         y = "Cell Count",
         fill = "Cell Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}



# subset seurat based on minimum number of cells filter:

subset_seurat_by_filtered_samples <- function(seurat_object, samples_to_keep) {
  if (length(samples_to_keep) > 0) {
    
    # Subset the Seurat object to keep only the samples that meet the criteria
    seurat_sub <- subset(seurat_object, subset = SampleID %in% samples_to_keep)
    
    return(seurat_sub)
  } else {
    print("No samples to keep.")
    return(NULL)
  }
}




## Violin plot with a list of samples to keep

plot_gene_expression <- function(seurat_object, samples_to_keep, gene, cell_types_of_interest) {
  if (length(samples_to_keep) > 0) {
    # Subset the Seurat object to keep only the samples that meet the criteria
    seurat_filtered <- subset(seurat_object, subset = SampleID %in% samples_to_keep)
    
    # Create a violin plot for the filtered data
    plot <- VlnPlot(seurat_filtered, features = gene, group.by = 'Class', idents = cell_types_of_interest, sort = TRUE, pt.size = 0.1) + ylim(c(0,5)) + NoLegend()
    
    return(plot)
  } else {
    print("No samples meet the criteria for plotting.")
  }
}



# filter_samples_by_gene_expression_O <- function(seurat_object, gene, threshold, percentage_threshold, cell_types_of_interest) {
#   samples_to_keep <- list()
#   sample_ids <- unique(seurat_object$SampleID)
#   seurat_object$SampleID <- as.character(seurat_object$SampleID)
#   
#   for (sample_id in sample_ids) {
#     results_per_cell_type <- list()
#     
#     for (cell_type in cell_types_of_interest) {
#       subset_cells <- subset(seurat_object, subset = SampleID == sample_id & Class == cell_type)
#       expression_data <- FetchData(subset_cells, vars = gene)
#       
#       if (nrow(expression_data) > 0) {  # Ensure there are cells to avoid division by zero
#         proportion_above_threshold <- sum(expression_data[, gene] > threshold, na.rm = TRUE) / nrow(expression_data) * 100
#         
#         if (proportion_above_threshold >= percentage_threshold) {
#           results_per_cell_type[[cell_type]] <- proportion_above_threshold
#         }
#       } else {
#         message(paste("No cells found for SampleID:", sample_id, "and Class:", cell_type))
#       }
#     }
#     
#     if (length(results_per_cell_type) > 0) {
#       samples_to_keep[[sample_id]] <- results_per_cell_type
#     }
#   }
#   
#   if (length(samples_to_keep) == 0) {
#     print("No samples met the threshold criteria.")
#   } else {
#     print(paste("Samples meeting criteria:", length(samples_to_keep)))
#   }
#   
#   return(names(samples_to_keep))
# }
# 
filter_samples_by_gene_expression_O <- function(seurat_object, gene, threshold, percentage_threshold, cell_types_of_interest) {
  samples_to_keep <- list()
  sample_ids <- unique(seurat_object$SampleID)
  seurat_object$SampleID <- as.character(seurat_object$SampleID)
  
  for (sample_id in sample_ids) {
    results_per_cell_type <- list()
    all_cell_types_met <- TRUE
    
    for (cell_type in cell_types_of_interest) {
      subset_cells <- subset(seurat_object, subset = SampleID == sample_id & Class == cell_type)
      expression_data <- FetchData(subset_cells, vars = gene)
      
      if (nrow(expression_data) > 0) {  # Ensure there are cells to avoid division by zero
        proportion_above_threshold <- sum(expression_data[, gene] > threshold, na.rm = TRUE) / nrow(expression_data) * 100
        
        if (proportion_above_threshold >= percentage_threshold) {
          results_per_cell_type[[cell_type]] <- list(
            proportion_above_threshold = proportion_above_threshold,
            num_cells = nrow(expression_data)
          )
        } else {
          all_cell_types_met <- FALSE
          break
        }
      } else {
        message(paste("No cells found for SampleID:", sample_id, "and Class:", cell_type))
        all_cell_types_met <- FALSE
        break
      }
    }
    
    if (all_cell_types_met) {
      samples_to_keep[[sample_id]] <- results_per_cell_type
    }
  }
  
  if (length(samples_to_keep) == 0) {
    print("No samples met the threshold criteria.")
  } else {
    print(paste("Samples meeting criteria:", length(samples_to_keep)))
  }
  
  for (sample_id in names(samples_to_keep)) {
    cat(paste("Sample ID:", sample_id, "\n"))
    for (cell_type in names(samples_to_keep[[sample_id]])) {
      cat(paste("  Cell Type:", cell_type, 
                "- Proportion Above Threshold:", samples_to_keep[[sample_id]][[cell_type]]$proportion_above_threshold,
                "- Number of Cells:", samples_to_keep[[sample_id]][[cell_type]]$num_cells, "\n"))
    }
  }
  
  return(names(samples_to_keep))
}

