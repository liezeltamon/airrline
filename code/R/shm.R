# From output of ddl.quantify_mutations(), calculate various mutation ratios:
# - CDR vs FWR mutation ratios
# - R vs S mutation ratios
# - By isotype (IGH, IGK, IGL)

calculate_mutation_ratios <- function(df) {
  
  # Helper function to safely calculate ratio
  safe_ratio <- function(numerator, denominator) {
    ratio <- numerator / denominator
    ratio[is.infinite(ratio) | is.nan(ratio)] <- NA
    return(ratio)
  }
  
  # Extract mutation type (count or freq) and region patterns
  mu_cols <- grep("^mu_(count|freq)_", names(df), value = TRUE)
  
  # Determine if we're working with counts or frequencies
  has_count <- any(grepl("^mu_count_", names(df)))
  has_freq <- any(grepl("^mu_freq_", names(df)))
  
  result_df <- df
  
  # Function to calculate ratios for a given metric type
  calc_ratios_for_type <- function(metric_type) {
    
    # Overall CDR vs FWR ratios
    # Total CDR (CDR1 + CDR2)
    cdr_total <- rowSums(df[, grep(paste0("^mu_", metric_type, "_cdr[12]_[rs]$"), 
                                    names(df), value = TRUE)], na.rm = TRUE)
    
    # Total FWR (FWR1 + FWR2 + FWR3)
    fwr_total <- rowSums(df[, grep(paste0("^mu_", metric_type, "_fwr[123]_[rs]$"), 
                                    names(df), value = TRUE)], na.rm = TRUE)
    
    result_df[[paste0("ratio_", metric_type, "_cdr_fwr")]] <<- safe_ratio(cdr_total, fwr_total)
    
    # R vs S mutation ratios
    # CDR R vs S
    cdr_r <- rowSums(df[, grep(paste0("^mu_", metric_type, "_cdr[12]_r$"), 
                                names(df), value = TRUE)], na.rm = TRUE)
    cdr_s <- rowSums(df[, grep(paste0("^mu_", metric_type, "_cdr[12]_s$"), 
                                names(df), value = TRUE)], na.rm = TRUE)
    result_df[[paste0("ratio_", metric_type, "_cdr_r_s")]] <<- safe_ratio(cdr_r, cdr_s)
    
    # FWR R vs S
    fwr_r <- rowSums(df[, grep(paste0("^mu_", metric_type, "_fwr[123]_r$"), 
                                names(df), value = TRUE)], na.rm = TRUE)
    fwr_s <- rowSums(df[, grep(paste0("^mu_", metric_type, "_fwr[123]_s$"), 
                                names(df), value = TRUE)], na.rm = TRUE)
    result_df[[paste0("ratio_", metric_type, "_fwr_r_s")]] <<- safe_ratio(fwr_r, fwr_s)
    
    # By isotype (IGH, IGK, IGL)
    for (isotype in c("IGH", "IGK", "IGL")) {
      # CDR vs FWR by isotype
      cdr_iso <- rowSums(df[, grep(paste0("^mu_", metric_type, "_cdr[12]_[rs]_", isotype, "$"), 
                                     names(df), value = TRUE)], na.rm = TRUE)
      fwr_iso <- rowSums(df[, grep(paste0("^mu_", metric_type, "_fwr[123]_[rs]_", isotype, "$"), 
                                     names(df), value = TRUE)], na.rm = TRUE)
      result_df[[paste0("ratio_", metric_type, "_cdr_fwr_", isotype)]] <<- safe_ratio(cdr_iso, fwr_iso)
      
      # R vs S in CDR by isotype
      cdr_r_iso <- rowSums(df[, grep(paste0("^mu_", metric_type, "_cdr[12]_r_", isotype, "$"), 
                                       names(df), value = TRUE)], na.rm = TRUE)
      cdr_s_iso <- rowSums(df[, grep(paste0("^mu_", metric_type, "_cdr[12]_s_", isotype, "$"), 
                                       names(df), value = TRUE)], na.rm = TRUE)
      result_df[[paste0("ratio_", metric_type, "_cdr_r_s_", isotype)]] <<- safe_ratio(cdr_r_iso, cdr_s_iso)
      
      # R vs S in FWR by isotype
      fwr_r_iso <- rowSums(df[, grep(paste0("^mu_", metric_type, "_fwr[123]_r_", isotype, "$"), 
                                       names(df), value = TRUE)], na.rm = TRUE)
      fwr_s_iso <- rowSums(df[, grep(paste0("^mu_", metric_type, "_fwr[123]_s_", isotype, "$"), 
                                       names(df), value = TRUE)], na.rm = TRUE)
      result_df[[paste0("ratio_", metric_type, "_fwr_r_s_", isotype)]] <<- safe_ratio(fwr_r_iso, fwr_s_iso)
    }
    
    # Individual region ratios
    # CDR1 vs FWR1, CDR2 vs FWR2, etc.
    for (region_num in 1:2) {
      cdr_region <- rowSums(df[, grep(paste0("^mu_", metric_type, "_cdr", region_num, "_[rs]$"), 
                                        names(df), value = TRUE)], na.rm = TRUE)
      fwr_region <- rowSums(df[, grep(paste0("^mu_", metric_type, "_fwr", region_num, "_[rs]$"), 
                                        names(df), value = TRUE)], na.rm = TRUE)
      result_df[[paste0("ratio_", metric_type, "_cdr", region_num, "_fwr", region_num)]] <<- 
        safe_ratio(cdr_region, fwr_region)
    }
  }
  
  # Calculate for counts if available
  if (has_count) {
    calc_ratios_for_type("count")
  }
  
  # Calculate for frequencies if available
  if (has_freq) {
    calc_ratios_for_type("freq")
  }
  
  return(result_df)
}

# Usage:
# df_with_ratios <- calculate_mutation_ratios(df)
