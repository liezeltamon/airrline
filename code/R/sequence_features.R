# Type of receptor sequence based on nucleotide and amino acid sequence

library(tidyverse)
library(alakazam)
library(immunarch)

# ----- FUNCTIONS -----

# df (input)  = Dataframe + name of column to use
# output = Dataframe of sample x features
# grouping_values = character vector to define grouping of repertoire sequence (e.g. sample id)

# 1. k-mer usage (nt or aa)

# immunarch::getKmers() counts k-mers using 1-nt sliding window as desired
tmp <- lapply(names(immdata$data), FUN = function(x) {
  cbind(grouping_id = rep(x), immdata$data[[x]])
})
df <- as.data.frame(do.call("rbind", tmp))

kmer_usage <- function(
    df, motif_length = 3, sequence_col, sequence_type, grouping_values, ...
) {
  # Needs both CDR3.nt and CDR3.aa columns because immunarch::getKmer() is checking internally
  # Create placeholder columns
  df["CDR3.nt"] <- "X"
  df["CDR3.aa"] <- "X"

  new_colname <- paste0("CDR3.", sequence_type)
  if (!identical(sequence_col, new_colname)) {
    df[[new_colname]] <- df[[sequence_col]]
    df[[sequence_col]] <- NULL
  }
  df_list <- split(df, f = grouping_values)
  #.coding Logical. If TRUE (by default) then removes all non-coding sequences from input data first.
  # Uses immunarch-specific columns
  kmer_counts_df <- immunarch::getKmers(
    df_list, motif_length, sequence_type, #...
  ) %>% 
    column_to_rownames("Kmer") %>% 
    t()
  kmer_fraction_df <- kmer_counts_df / rowSums(kmer_counts_df)
  kmer_fraction_df
}
output_df <- kmer_usage(
  df, motif_length = 3, sequence_col = "CDR3.nt", sequence_type = "nt", grouping_values = df$grouping_id
)

# 2. aa motif usage: GLIPH2 + python processing

# 3. Physiochemical properties

#df <- as.data.frame(ExampleDb)
output_df <- alakazam::aminoAcidProperties(
  df, seq = "junction", nt = TRUE, trim = FALSE, label = "amino_acid_properties"
) %>% 
  select(contains("amino_acid_properties"))
