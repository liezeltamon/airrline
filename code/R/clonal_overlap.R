# Receptor sequence overlap between cell/sample subsets (e.g. cell types, isotypes, tissues)

# Setup

library(assertthat)
library(tidyverse)
library(scRepertoire)

# ----- FUNCTIONS -----

.clonal_overlap <- function(
    x,                    # df with clone_key and subset_key columns
    subset_values,        # Subset values of cells in x in same order
    subset_levels = NULL, # All possible subset levels
    method = c("overlap", "morisita", "jaccard", "cosine", "raw"), # As in scRepertoire:::clonalOverlap()
    normalise = FALSE     # Normalise final values? Depends on the overlap metric and may not be necessary after subsampling
) {
  
  # Checks
  method <- match.arg(method)
  assert_that(
    length(x) == length(subset_values),
    msg = ".clonal_overlap(): Length of x and subset_values not equal"
  )
  
  # scRepertoire:::.calculateIndex expects a list of dataframe with clone ids as a column
  # Variables, e.g. cloneCall, input.data, length, used as in scRepertoire:::clonalOverlap(), 
  # whose internal functions are being used here
  cloneCall = "clone_id"
  x_df <- data.frame(df$clone_id)
  colnames(x_df) <- cloneCall
  if(method == "morisita") {
    return_type <- "freq"
  } else {
    return_type <- "unique"
  }
  
  input.data <- split(x_df, f = subset_values)
  num_samples <- length(input.data)
  names_samples <- names(input.data)
  length <- seq_len(num_samples)
  
  # Choose overlap metric function
  indexFunc <- switch(
    method,
    "morisita" = scRepertoire:::.morisitaCalc,
    "jaccard"  = scRepertoire:::.jaccardCalc,
    "raw"      = scRepertoire:::.rawCalc,
    "overlap"  = scRepertoire:::.overlapCalc,
    "cosine"  =  scRepertoire:::.cosineCalc,
    stop("Invalid method provided")
  )
  
  coef_matrix <- data.frame(matrix(NA, num_samples, num_samples))
  coef_matrix <- scRepertoire:::.calculateIndex(
    input.data,
    length,
    cloneCall,
    coef_matrix, 
    indexFunc, 
    return_type
  )
  colnames(coef_matrix) <- names_samples
  rownames(coef_matrix) <- names_samples
  
  output <- coef_matrix %>% 
    rownames_to_column("row_id") %>%
    pivot_longer(-row_id, names_to = "col_id") %>% 
    filter(!is.na(value)) %>% 
    unite("features", row_id, col_id, sep = ".") %>% 
    pivot_wider(names_from = features, values_from = value) %>% 
    # 1-row dataframe to named numeric vector to match .inter_subset output
    as.list() %>% unlist()
  
  if (normalise) {
    output <- (output - min(output)) / (max(output) - min(output))
  }
  
  return(output)
}
