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

.grouping_size <- function(
  x,
  ...
) {
  return(length(x))
}
