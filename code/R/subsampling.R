# Re-implementation of scClonetoire by Rachael Bashford-Rogers, Sakina Amin, et al.
# Made to work with a per-cell data frame as input, intergrated inter- and intra in one function
# added parallelisation with foreach, seeding for reproducibility
# Other utility functions to help choose subsample size, classify clones based on size
# Other clonality measures needing subsampling can easily be integrated as additional methods

library(assertthat)
library(foreach)
library(doParallel)
library(tidyverse)
#source("misc.R")
#source("subset_clonality.R")

#source("utils/env.R")
#library(scRepertoire) # For some diversity measures

# Plot distribution of grouping sizes to determine ideal size of a boot/subsample
get_grouping_sizes <- function(
    data,
    grouping_keys,
    export_table = FALSE, # Match screpertoire argument
    drop_na = FALSE, # Don't count rows in data with NA in any columns in grouping_keys
    bins = 30
) {
  
  p_df <- data %>% 
    group_by(across(all_of(grouping_keys))) %>% 
    count() %>% 
    arrange(n) %>% 
    { if (drop_na) drop_na(.) else . }
  
  if (export_table) {
    return(p_df)
  } else {
    warning("Apply heuristic to choose.")
    p <- ggplot(p_df, aes(n)) +
      geom_histogram(bins = bins) +
      theme_bw()
    print(p)
  }
  
}

measure_clonality <- function(
    data, # Per-cell dataframe, separate B/CD4/CD8T so can also submit as separate jobs
    grouping_keys, # Columns to determine group to be subsampled
    subset_key = NULL, # For inter-subset clonality
    clone_key,
    drop_na = FALSE,
    method_fun = NULL, # Ideally written like other pre-defined method functions giving 1 value result
    method = "intra_subset", # Could be diversity metrics or any other metric that can be measured per boot/subsample
    # Set NULL if not subsampling; grouping not considered for calculation if size below this
    subsample_size = NULL,
    n_subsamples = 1,
    seed = 42,
    n_cpu = 1, # Parallelisation is per grouping
    return.boots = TRUE, # As in scRepertoire, return all calculations otherwise summarise subsampled values 
    sep = "_", # When merging grouping_keys to grouping_id column
    ... # Additional arguments to pre-defined method functions
) {

  # Assign right function to calculate_measure() based on metric argument, use switch?
  # See dedicated scripts for functions to assign to .calculate_measure() applying to vector of clone ids (x)
  if (!is.null(method_fun)) {
    .calculate_measure <- method_fun
    method <- NULL
  } else {
    .calculate_measure <- switch(
      method,
      "clone_size"     = .clone_size,
      "grouping_size"  = .grouping_size,
      "intra_subset"   = .intra_subset,
      "inter_subset"   = .inter_subset,
      # Add diversity metrics here
      stop("measure_clonality(): Invalid method provided")
    )
  }

  keys_to_select <- c("grouping_id", clone_key, subset_key)
  keys_to_select <- keys_to_select[!is.null(keys_to_select)]
  data_wrangled <- data %>%
    unite("grouping_id", grouping_keys, remove = FALSE, sep = sep) %>% 
    #select(all_of(c("grouping_id", clone_key, subset_key))) %>% 
    select(all_of(keys_to_select)) %>%
    rename_with(~ "clone_id", all_of(clone_key)) %>% 
    { if (drop_na) drop_na(.) else . }
  
  # List of df per group
  grouped_data_lst <- split(
    data_wrangled,
    data_wrangled[ ,"grouping_id", drop = FALSE]
  )
  
  # Get required information across groupings

  # -Get seed values for each group to ensure reproducibility across runs and parallelisation
  set.seed(seed)
  seed_values <- sample(
    1:(.Machine$integer.max), size = length(grouped_data_lst), replace = FALSE
  )
  
  # -Subsample size for each group, either set by user or based on group size (if subsample_size = NULL)
  if (!is.null(subsample_size)) {
    subsample_sizes <- rep(subsample_size, times = length(grouped_data_lst))
  } else {
    subsample_sizes <- sapply(grouped_data_lst, nrow)
  }
  rm(subsample_size)
  
  # -Subset levels for each group, mainly used for inter-subset clonality
  subset_levels <- "value"
  if (!is.null(subset_key)) {
    subset_levels <- sort(unique(data[[subset_key]]))
  }
  
  # -Add these to grouped_data_lst so they are passed to workers when parallelising, also add group name for easier debugging
  for (i in seq_along(grouped_data_lst)) {
    grouped_data_lst[[i]] <- list(
      df = grouped_data_lst[[i]],
      seed = seed_values[i],
      name = names(grouped_data_lst)[[i]],
      subset_levels = subset_levels,
      subsample_size = subsample_sizes[i]
    )
  }

  ##### parallel start
  
  if(n_cpu > 1){
    #registerDoParallel(cores = 2) # For some reason, outputs collected from workers wrong when using this
    cl <- makeCluster(n_cpu)
    registerDoParallel(cl)
    `%op%` <- `%dopar%`
    print(paste0("Running with ", n_cpu, " cores."), quote = FALSE)
  } else {
    `%op%` <- `%do%`
  }
  
  SUBSAMPLED <- foreach(
    grp_data = grouped_data_lst,
    .inorder = TRUE,
    .packages = c("tidyverse", "assertthat"),
    .export = c("n_subsamples", "subset_key", ".calculate_measure")
  ) %op% {
    
    message(
      grp_data$name, " subsample_size=", grp_data$subsample_size, " seed=", grp_data$seed
    )
    
    if (nrow(grp_data$df) >= grp_data$subsample_size) {
      subsampled_lst <- vector("list", length = n_subsamples)
      set.seed(grp_data$seed)
      for (b in 1:n_subsamples) {
        sample_inds <- sample(
          1:nrow(grp_data$df), size = grp_data$subsample_size, replace = FALSE
        )
        input <- grp_data$df[sample_inds, "clone_id"]
        subset_values <- if (!is.null(subset_key)) grp_data$df[sample_inds, subset_key] else NULL
        # subset_values ignored for other functions except .inter_subset()
        subsampled_lst[[b]] <- .calculate_measure(
          input,
          subset_values = subset_values,
          subset_levels = grp_data$subset_levels
        ) %>% 
        # Convert to data.frame if not already, to match output of .inter_subset() and make downstream summarisation easier
        { if (is.data.frame(.)) . else data.frame(value = .) }
      }
      message(grp_data$name, " subsampling done!")
    } else {
      subsampled_lst <- list(as.data.frame(
        matrix(NA, nrow = 1, ncol = length(grp_data$subset_levels), dimnames = list(NULL, grp_data$subset_levels))
      ))
      message(
        grp_data$name, " has size ", nrow(grp_data$df), " < subsample_size ", grp_data$subsample_size, 
        " so returning NA."
      )
    }
    
    message(grp_data$name, " done!")
    
    return(subsampled_lst)
    
  }
  
  if(n_cpu > 1){
    stopCluster(cl)
  }
  
  ##### parallel end
  
  assert_that(length(grouped_data_lst) == length(SUBSAMPLED))
  names(SUBSAMPLED) <- names(grouped_data_lst)
  
  # Return or summarise subsampled values (one value per grouping)
  
  if (return.boots) {
    OUTPUT <- .transform_to_df(SUBSAMPLED, method = method)
  } else {
    .summarise_measure <- .summarise_one_value
    if (!is.null(method) && method == "inter_subset") {
        .summarise_measure <- .summarise_inter_subset
    }
    OUTPUT <- cbind.data.frame(
      grouping_id = names(grouped_data_lst),
      value = .summarise_measure(SUBSAMPLED)
    )
  }

  return(OUTPUT)
}

# Functions to wrangle / summarise to get OUTPUT:

# Convert SUBSAMPLED to melted data.frame ready for plotting
# grouping_id with NA as value as those repertoire with size < subsample_size
.transform_to_df <- function(x, method = NULL) {
  grp_df_lst <- sapply(x, USE.NAMES = TRUE, simplify = FALSE, function(df_lst) {
    bind_rows(df_lst)
  })
  df <- bind_rows(grp_df_lst, .id = "grouping_id")
  if (!is.null(method) && method == "inter_subset") {
    df <- df %>% pivot_longer(!grouping_id, names_to = "subset_id")
  }
  return(df)
}

# For measures that return one value per subsample
.summarise_one_value <- function(lst) {
  # One final value per grouping
  unlist(
    lapply(lst, function(x) mean(unlist(x, recursive = TRUE)))
  )
}

.summarise_inter_subset <- function(lst) {
  # One final value per grouping, subset_id combination
  bind_rows(
    y <- lapply(lst, function(x) {
      x_df <- bind_rows(x)
      colMeans(x_df)
    })
  )
}

x <- SUBSAMPLED[[1]]
with_value <- lapply(y, function(x){
  if ("value" %in% names(x)) {
    return(TRUE)
  } else {
    return(FALSE)
  } 
}) %>% unlist()