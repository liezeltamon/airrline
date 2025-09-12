# VDJ gene usage: Count and relative abundance of V(D)J alleles, genes or families within groups

library(tidyverse)
library(alakazam)

# ----- MAIN -----

#df <- ExampleDb
df <- as.data.frame(vdj)
df$grouping_id <- "a"
output_df <- alakazam::countGenes(
  df, gene = "v_call_genotyped", groups = "grouping_id", mode = "allele"
)
output_df <- output_df %>% 
  pivot_wider(
    id_cols = "grouping_id", names_from = "gene", values_from = "seq_freq"
  ) %>% 
  column_to_rownames("grouping_id")
