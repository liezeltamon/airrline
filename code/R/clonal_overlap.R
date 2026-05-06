# Receptor sequence overlap between cell/sample subsets (e.g. cell types, isotypes, tissues)


# set.seed(42)

# # Simulate clone ids - some shared between subsets to make overlap meaningful
# n_cells <- 200
# clone_pool <- paste0("clone_", 1:30)  # 30 unique clones

# x <- sample(clone_pool, n_cells, replace = TRUE, prob = c(
#   rep(0.1, 5),   # 5 dominant clones
#   rep(0.02, 25)  # 25 rare clones
# ))

# # Subset values - 3 cell types with uneven sizes
# subset_values <- sample(c("t.cd4", "t.cd8", "t.reg"), n_cells, replace = TRUE, prob = c(0.5, 0.35, 0.15))


# Setup

library(assertthat)
library(tidyverse)
library(scRepertoire)

# ----- FUNCTIONS -----

.clonal_overlap <- function(
    x,              # Clone ids
    subset_values,  # Subset values of cells in x in same order
    subset_levels,  # all possible subset levels
    subset_levels_keep = NULL, # If not NULL, only keep overlap values for the specified subset pairs, set rest to NA.
    # e.g. 
    # list(
    #   c("t1__naive", "t2__treg"),
    #   c("t1__naive", "t2__memory"),
    #   c("t1__memory", "t2__treg")
    # )
    method = c("raw", "overlap", "morisita", "jaccard", "cosine"), # As in scRepertoire:::clonalOverlap(). - see https://github.com/BorchLab/scRepertoire/blob/HEAD/R/clonalOverlap.R
    fill_value = 0, # Value used for missing subset pairs; meaning depends on the overlap metric
    normalise = FALSE, # Will only apply for now if method == "raw"
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
    msg = paste0(
      ".clonal_overlap(): Not all subset_values are in subset_levels (supposed to be the subset universe). ",
      "Check subset_values not in subset_levels: ",
      paste(setdiff(unique(subset_values), subset_levels), collapse = ", ")
    )
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
    fill_value, # Use a non-NA fill so all possible subset pairs based on subset_levels remain represented in the output.
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

  coef_matrix <- coef_matrix[
    sort(rownames(coef_matrix), decreasing = FALSE),
    sort(colnames(coef_matrix), decreasing = FALSE)
  ]

  output <- coef_matrix %>% 
    as.data.frame() %>%
    rownames_to_column("row_id") %>%
    pivot_longer(-row_id, names_to = "col_id") %>% 
    filter(!is.na(value)) %>% 
    
    # Standardise subset-pair names by sorting levels alphabetically first.
    rowwise() %>%
    mutate(
      features = paste(sort(c(row_id, col_id)), collapse = ".."),
      row_id = NULL, col_id = NULL
    ) %>%
    ungroup() %>%
    #unite("features", row_id, col_id, sep = "..") %>% 

    pivot_wider(names_from = features, values_from = value) %>% 
    # 1-row dataframe to named numeric vector to match .inter_subset output
    as.list() %>%
    unlist()
  
  # Keep only the requested subset pairs and set all other overlap features to NA.
  if (!is.null(subset_levels_keep)) {
    subset_levels_keep_colnames <- vapply(
      subset_levels_keep,
      function(pair) {
        paste(sort(pair), collapse = "..")
      },
      character(1)
    )

    missing_pairs <- subset_levels_keep_colnames[
      !subset_levels_keep_colnames %in% names(output)
    ]
    if (length(missing_pairs) > 0) {
      stop(
        paste0(
          ".clonal_overlap(): The following subset pairs in ",
          "subset_levels_keep are missing from output: ",
          paste(missing_pairs, collapse = ", ")
        )
      )
    }

    output[!(names(output) %in% subset_levels_keep_colnames)] <- NA_real_
    message(
      ".clonal_overlap(): Returning overlap only for subset pairs in subset_levels_keep; other pairs set to NA."
    )
  }

  if (normalise) {
    if (method == "raw") {
      output_non_na <- output[!is.na(output)]
      stopifnot(length(output_non_na) > 0)
      # Normalise to sum to 1
      if (sum(output_non_na) > 0) {
        output[!is.na(output)] <- output_non_na / sum(output_non_na)
      }
    } else {
      warning("Normalisation can only be applied to raw method for now")
    }
  }
  
  n_expected_pairs <- choose(length(subset_levels), 2)
  if (length(output) != n_expected_pairs) {
    stop(
      ".clonal_overlap(): output length ", length(output),
      " != expected pair count ", n_expected_pairs
    )
  }

  return(output)
}
