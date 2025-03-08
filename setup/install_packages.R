version # Print R version
renv::load()
mass_version = "7.3-60"
matrix_version = "1.6-5"

#

renv::install(
  
  c(
    "BiocManager",
    "rprojroot",
    "devtools",
    paste0("MASS@", mass_version),
    paste0("Matrix@", matrix_version),
    "tidyverse",
    "doParallel",
    "foreach",
    "cowplot",
    "pheatmap",
    "svglite",   # Loaded by utils/env.R for saving svg ggplots
    "reticulate" # Loaded by utils/env.R to use conda/mamba envs but consider removing
  )
  
)
