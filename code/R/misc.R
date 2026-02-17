# Functions for other repertoire measures

.clone_size <- function(
    x,
    scale = "auto",
    ...
) {
  tbl <- table(x)
  tbl_numeric <- as.numeric(tbl)
  names(tbl_numeric) <- names(tbl)
  if (scale == "auto") {
    output <- tbl_numeric / length(x)
  } else if (scale == "raw") {
    output <- tbl_numeric
  } else if (is.numeric(scale)) {
    output <- tbl_numeric / scale
  } else {
    stop(".clone_size(): Invalid scale argument. Indicate 'auto', 'raw' or a numeric value.")
  }
  return(as.data.frame(enframe(output, name = "clone_id", value = "value")))
}

.grouping_size <- function(x, ...) {
  return(length(x))
}

# Repertoire convergence (currently applicable for single-cell,
# For bulk need to deal with read signal or normalised reads)
# Needs subsampling because ratio is sensitive to raw number of unique sequences
# **Need to adapt subsampling function for this**
# The lower the value the more convergent
# Repertoire being convergent actually requires more resources so a repertoire being like that means that it is doing it for a purpose
repertoire_convergence <- function (x_nt, x_aa) {
  n_unique_nt_seq <- length(unique(x_nt))
  n_unique_aa_seq <- length(unique(x_aa))
  n_unique_aa_seq / n_unique_nt_seq
}
