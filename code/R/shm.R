# From output of ddl.quantify_mutations(), calculate additional mutation features

.compute_ratio <- function(numerator, denominator, pseudocount = 1) {
  ratio <- log2((numerator + pseudocount) / (denominator + pseudocount))
  n_nonfinite <- sum(!is.finite(ratio))
  stopifnot("Non-finite values in ratio computation" = n_nonfinite == 0)
  ratio
}

.make_name_ratio <- function(num, den, prefix) paste(
  prefix, num, "over", den, sep = "_"
)

.compute_ratio_combinations <- function(
  numerator_mx, denominator_mx, prefix = "muratio", pairwise = TRUE
) {
  stopifnot(is.matrix(numerator_mx), is.matrix(denominator_mx))
  num_names <- colnames(numerator_mx)
  den_names <- colnames(denominator_mx)
  nrow_out <- nrow(numerator_mx)

  if (pairwise) {
    combos <- expand.grid(num = num_names, den = den_names, stringsAsFactors = FALSE)
  } else {
    stopifnot(length(num_names) == length(den_names))
    combos <- data.frame(num = num_names, den = den_names, stringsAsFactors = FALSE)
  }

  ncol_out <- nrow(combos)
  output <- matrix(NA_real_, nrow = nrow_out, ncol = ncol_out)
  out_names <- character(ncol_out)

  for (i in seq_len(ncol_out)) {
    num <- combos$num[i]
    den <- combos$den[i]
    output[, i] <- .compute_ratio(numerator_mx[, num], denominator_mx[, den])
    out_names[i] <- .make_name_ratio(num, den, prefix = prefix)
  }
  colnames(output) <- out_names
  output
}

.parse_colnames <- function(colnames_vec) {
  # Splits column names by "_" and returns a data.frame with each part as a column
  # Handles missing parts by filling with NA
  split_parts <- strsplit(colnames_vec, "_")
  max_len <- max(lengths(split_parts))
  # Pad each split to max_len with NA
  padded <- lapply(split_parts, function(x) { length(x) <- max_len; x })
  parts_df <- as.data.frame(do.call(rbind, padded), stringsAsFactors = FALSE)
  colnames(parts_df) <- paste0("part", seq_len(max_len))
  rownames(parts_df) <- colnames_vec
  # Change colnames
  colnames(parts_df) <- c("prefix", "value_type", "region", "mutation_type", "chain")[1:max_len]

  # Split part3 (e.g., fwr1, cdr2) into region and number
  if ("region" %in% colnames(parts_df)) {
    parts_df$region_type <- gsub("[0-9]+$", "", parts_df$region)
    parts_df$region_n <- gsub("^[a-zA-Z]+", "", parts_df$region)
  }
  
  parts_df
}

compute_shm_features <- function(df) {
  
  feature_list <- list()
  parse_df <- .parse_colnames(names(df))

  # Groupings per value_type, mutation_type, chain

  groupings_df <- na.omit(unique(
    parse_df[, c("value_type", "mutation_type", "chain")]
  ))

  for (g in seq_len(nrow(groupings_df))) {

    grouping <- groupings_df[g, ]

    # 1. Calculate pairwise fwr/cdr ratios

    cdr_cols_g <- rownames(
      subset(
        parse_df,
        value_type == grouping$value_type &
        mutation_type == grouping$mutation_type &
        chain == grouping$chain &
        region_type == "cdr"
      )
    ) %>% sort()

    fwr_cols_g <- rownames(
      subset(
        parse_df,
        value_type == grouping$value_type &
        mutation_type == grouping$mutation_type &
        chain == grouping$chain &
        region_type == "fwr"
      )
    ) %>% sort()

    mx <- .compute_ratio_combinations(
      df[ , fwr_cols_g, drop = FALSE] %>% as.matrix(),
      df[ , cdr_cols_g, drop = FALSE] %>% as.matrix(),
      prefix = "muratio"
    )
    rownames(mx) <- rownames(df)
    feature_list <- c(feature_list, list(mx))
  
    # 2. Average fwr regions

    fwr_ave <- apply(
      df[ , fwr_cols_g, drop = FALSE], 1, mean, na.rm = FALSE
    )

    # 3. Average fwr and cdr regions

    .compute_average <- function(df, cols) {
      feature_name <- paste0(
        "mumean", "_", unique(gsub("[0-9]+", "", cols))
      )
      stopifnot(length(feature_name) == 1)
      setNames(
        list(apply(df[, cols, drop = FALSE], 1, mean, na.rm = FALSE)),
        feature_name
      )
    }
    cdr_ave <- .compute_average(df, cdr_cols_g)
    fwr_ave <- .compute_average(df, fwr_cols_g)
    feature_list[[names(cdr_ave)]] <- cdr_ave[[1]]
    feature_list[[names(fwr_ave)]] <- fwr_ave[[1]]
    
    # 4. Average fwr regions / average cdr regions

    feature_name <- .make_name_ratio(
      names(fwr_ave), names(cdr_ave), prefix = "muratio"
    )
    feature_list[[feature_name]] <- .compute_ratio(
      fwr_ave[[1]], cdr_ave[[1]]
    )
  }

  # Groupings per value_type, chain, region_type, region_n

  groupings_df <- na.omit(unique(
    parse_df[, c("value_type", "chain", "region_type", "region_n")]
  ))

  for (g in seq_len(nrow(groupings_df))) {

    grouping <- groupings_df[g, ]

    # 5. Ratio s / r in each region

    cols_g <- rownames(
      subset(
        parse_df,
        value_type == grouping$value_type &
        chain == grouping$chain &
        region_type == grouping$region_type &
        region_n == grouping$region_n
      )
    ) %>% sort()
    stopifnot(length(cols_g) == 2)
    feature_name <- .make_name_ratio(
      cols_g[grepl("_r_", cols_g)],
      cols_g[grepl("_s_", cols_g)],
      prefix = "muratio"
    )
    feature_list[[feature_name]] <- .compute_ratio(
      df[, cols_g[grepl("_r_", cols_g)]],
      df[, cols_g[grepl("_s_", cols_g)]]
    )
  }

  do.call(cbind, feature_list)

}

# Example:
#   df <- data.frame(...)
#   compute_shm_features(df)
