# Measure clonality of cell subsets e.g. cell types, isotypes

library(tidyverse)
library(assertthat)

# ----- FUNCTIONS -----

.intra_subset <- function(
    x,                     # a subsample of clone ids
    expanded_min_size = 2, # clone size
    scale = "auto",
    ...                    # to ignore .inter_subset() arguments passed to sure()
) { 
  tbl <- table(x)
  n_in_expanded <- sum(tbl[tbl >= expanded_min_size])
  if (scale == "auto") {
    output <- n_in_expanded / length(x)
  } else if (scale == "raw") {
    output <- n_in_expanded
  } else if (is.numeric(scale)) {
    output <- n_in_expanded / scale
  } else {
    stop(".intra_subset(): Invalid scale argument. Indicate 'auto', 'raw' or a numeric value.")
  }
  return(output)
}

.inter_subset <- function(
    x,                    # df with clone_key and subset_key columns
    subset_values,        # subset values of cells in x in same order
    subset_levels = NULL, # all possible subset levels
    method_inter = c("sum_n", "average_proportions"),
    expanded_min_size = 3 # min size of clone to be considered expanded
) {
  
  # Prepare functions
  .sum_n <- function(tbl) {
    # Original implementation does scaling using all subsamples i.e.:
    # colSums(subsamples) / sum(subsamples)
    # But i think below version account for variation across subsamples which hopefully should be minimal
    return(colSums(tbl) / sum(tbl))
  }
  .average_proportions <- function(tbl) {
    prop_tbl <- apply(tbl, MARGIN = 1, proportions)
    return(apply(prop_tbl, MARGIN = 2, mean))
  }
  
  method_inter <- method_inter[1]
  
  assert_that(
    length(x) == length(subset_values),
    msg = ".inter_subset(): Length of x and subset_values not equal"
  )
  
  .calculate_output <- switch(
    method_inter,
    "sum_n" = .sum_n,
    "average_proportions" = .average_proportions,
    stop(".inter_subset(): Invalid method provided")
  )
  
  # Converted to factor so all levels are used
  x_fct <- factor(x)
  subset_values_fct <- factor(subset_values, levels = subset_levels)
  tbl <- table(x_fct, subset_values_fct)
  expanded_clones <- which(rowSums(tbl) >= expanded_min_size)
  
  if (length(expanded_clones) >= 1) {
    output <- .calculate_output(tbl[expanded_clones, , drop = FALSE])
    # To match format below
    output <- t(output) %>% as.data.frame()
  } else {
    output <- as.data.frame(
      matrix(data = 0, nrow = 1, ncol = ncol(tbl),
             dimnames = list("", colnames(tbl)))
    )
  }
  return(output)
}
