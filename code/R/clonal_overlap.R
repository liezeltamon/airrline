# Receptor sequence overlap between cell/sample subsets (e.g. cell types, isotypes, tissues)


set.seed(42)

# Simulate clone ids - some shared between subsets to make overlap meaningful
n_cells <- 200
clone_pool <- paste0("clone_", 1:30)  # 30 unique clones

x <- sample(clone_pool, n_cells, replace = TRUE, prob = c(
  rep(0.1, 5),   # 5 dominant clones
  rep(0.02, 25)  # 25 rare clones
))

# Subset values - 3 cell types with uneven sizes
subset_values <- sample(c("t.cd4", "t.cd8", "t.reg"), n_cells, replace = TRUE, prob = c(0.5, 0.35, 0.15))


# Setup

library(assertthat)
library(tidyverse)
library(scRepertoire)

# ----- FUNCTIONS -----

.clonal_overlap <- function(
    x,              # Clone ids
    subset_values,  # Subset values of cells in x in same order
    subset_levels,  # all possible subset levels
    method = c("overlap", "morisita", "jaccard", "cosine", "raw"), # As in scRepertoire:::clonalOverlap()
    normalise = FALSE,     # Normalise final values? Depends on the overlap metric and may not be necessary after subsampling
    ...
) {
  
  # Checks
  method <- match.arg(method)
  assert_that(
    length(x) == length(subset_values),
    msg = ".clonal_overlap(): Length of x and subset_values not equal"
  )
  assert_that(
    all(subset_values %in% subset_levels),
    msg = ".clonal_overlap(): Not all subset_values are in subset_levels (supposed to be the subset universe)"
  )
  assert_that(
    length(unique(subset_values)) >= 2,
    msg = ".clonal_overlap(): Need at least 2 subsets to calculate overlap"
  )

  # scRepertoire:::.calculateIndex expects a list of dataframe with clone ids as a column
  # Variables, e.g. cloneCall, input.data, length, used as in scRepertoire:::clonalOverlap(), 
  # whose internal functions are being used here
  cloneCall = "clone_id"
  x_df <- data.frame(x)
  colnames(x_df) <- cloneCall
  if(method == "morisita") {
    return_type <- "freq"
  } else {
    return_type <- "unique"
  }
  
  input.data <- split(x_df, f = subset_values)
  #input.data <- input.data[sort(names(input.data), decreasing = FALSE)]
  #subset_levels <- sort(subset_levels, decreasing = FALSE)
  num_samples <- length(input.data)
  names_samples <- names(input.data)
  idx_samples <- seq_len(num_samples)
  
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

  # Initialise empty matrix with all possible pairs based on subset_levels, to be filled in by scRepertoire:::.calculateIndex() and returned as output
  subset_levels <- union(names(input.data), subset_levels)
  num_subset_levels <- length(subset_levels)
  coef_matrix <- matrix(
    0, # 0 vs. NA so all subset levels considered as pairs with 0 overlap instead of being ignored, then all pairs not needed to be returned set to NA later (non-lower tri)
    num_subset_levels, num_subset_levels,
    dimnames = list(
      # Cause scRepertoire:::.calculateIndex will only fill levels = names(input.data) in that order 
      subset_levels, subset_levels
    )
  )

  coef_matrix <- scRepertoire:::.calculateIndex(
    input.data,
    idx_samples,
    cloneCall,
    coef_matrix,
    indexFunc,
    return_type
  )
  coef_matrix[lower.tri(coef_matrix, diag = TRUE)] <- NA_real_

  # To enforce that pair names are formed by sorting alphabetically first, e.g. A_B not B_A, making names of output predictable and can be created outside the function
  coef_matrix <- coef_matrix[
    sort(rownames(coef_matrix), decreasing = FALSE),
    sort(colnames(coef_matrix), decreasing = FALSE)
  ]

  output <- coef_matrix %>% 
    as.data.frame() %>%
    rownames_to_column("row_id") %>%
    pivot_longer(-row_id, names_to = "col_id") %>% 
    filter(!is.na(value)) %>% 
    unite("features", row_id, col_id, sep = "..") %>% 
    pivot_wider(names_from = features, values_from = value) %>% 
    # 1-row dataframe to named numeric vector to match .inter_subset output
    as.list() %>% unlist()
  
  if (normalise) {
    rng <- max(output) - min(output)
    if (rng > 0) output <- (output - min(output)) / rng
    else stop(".clonal_overlap(): All overlap values identical; returning unnormalised values.")
  }
  
  return(output)
}
