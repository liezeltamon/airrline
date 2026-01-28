library(immApex)
library(assertthat)
library(Biostrings)

is_valid_sequence <- function(sequences, type = c("nt", "aa"), productive_only = FALSE) {
  type <- match.arg(type)
  
  sapply(sequences, function(seq) {
    if (is.na(seq) || seq == "" || !is.character(seq)) {
      return(FALSE)
    }
    
    # Check for stop codons if productive_only is TRUE
    if (type == "aa" && productive_only && grepl("\\*", seq)) {
      return(FALSE)
    }
    
    # Validate with Biostrings
    tryCatch({
      if (type == "nt") {
        Biostrings::DNAString(seq)
      } else if (type == "aa") {
        Biostrings::AAString(seq)
      }
      return(TRUE)
    }, error = function(e) {
      return(FALSE)
    })
  })
}

get_all_kmers <- function(k, type = c("nt", "aa")) {
  type <- match.arg(type)
  
  if (type == "nt") {
    as.character(Biostrings::mkAllStrings(DNA_BASES, k))
  } else if (type == "aa") {
    as.character(Biostrings::mkAllStrings(AA_STANDARD, k))
  } else {
    stop("Invalid type. Choose 'nt' or 'aa'.")
  }
}

kmer_usage <- function(
    df,
    sequence_key,
    grouping_key = NULL,
    motif_length = 3,
    sequence_type = c("nt", "aa"),
    output_type = c("percentage", "fraction", "count"),
    min.depth = 1,
    ...
) {
  
  warning("This function cannot deal with invalid sequences so deal with them beforehand.")

  # Validate inputs
  assertthat::assert_that(is.data.frame(df))
  assertthat::assert_that(sequence_key %in% names(df))
  # Assert that sequences are productive i.e. not containing "*"
  assertthat::assert_that(
    all(!grepl("\\*", df[[sequence_key]])),
    msg = "Sequences contain stop codons ('*'). Please filter them out before proceeding."
  )
  if (!is.null(grouping_key)) {
    assertthat::assert_that(grouping_key %in% names(df))
  }
  assertthat::assert_that(motif_length > 0)

  # Match arguments
  sequence_type <- match.arg(sequence_type)
  output_type <- match.arg(output_type)
  
  # Prepare input
  if (!is.null(grouping_key)) {
    # Make list of sequences per group
    input <- split(df[[sequence_key]], df[[grouping_key]])
    output_row_ids <- names(input)
  } else {
    input <- df[[sequence_key]]
    # Sequence can have duplicate entries, so use row numbers as IDs
    output_row_ids <- seq_len(nrow(df))
  }

  # Count k-mers per sequence/sample
  kmer_counts_list <- sapply(
    input,
    USE.NAMES = TRUE,
    simplify = FALSE,
    FUN = immApex:::calculateMotif,
    motif.lengths = motif_length,
    min.depth = min.depth,
    ...
  )
  
  # Convert kmer_counts_list to sequence/sample x all possible k-mers matrix
  
  # Vector of all possible k-mers of given length
  all_kmers <- get_all_kmers(motif_length, sequence_type)

  kmer_mx <- matrix(
    0,
    nrow = length(kmer_counts_list),
    ncol = length(all_kmers),
    dimnames = list(seq_along(kmer_counts_list), all_kmers)
  )
  
  # Fill in counts
  for (i in seq_along(kmer_counts_list)) {
    counts_df <- kmer_counts_list[[i]]
    kmer_mx[i, counts_df$motif] <- counts_df$frequency
    message("Processed ", i, " / ", length(kmer_counts_list), " entries.")
  }

  # Convert to requested output format
  if (output_type == "percentage") {
    kmer_mx <- kmer_mx / rowSums(kmer_mx) * 100
    assertthat::assert_that(all(abs(rowSums(kmer_mx) - 100) < 1e-10))
  } else if (output_type == "fraction") {
    kmer_mx <- kmer_mx / rowSums(kmer_mx)
    assertthat::assert_that(all(abs(rowSums(kmer_mx) - 1) < 1e-10))
  } else if (output_type == "count") {
    # Do nothing, already counts
  } else {
    stop("Invalid output_type. Choose 'count', 'percentage', or 'fraction'.")
  }

  assertthat::assert_that(nrow(kmer_mx) == length(input))
  rownames(kmer_mx) <- output_row_ids
  kmer_mx
}
# #Example usage:
# df <- data.frame(
#   sequence = c("ATCGATCG", "GCTAGCTA", "ATCGATCG", "GCTAGCTA"),
#   group = c("A", "A", "B", "B")
# )
# result <- kmer_usage(
#   df,
#   sequence_key = "sequence",
#   grouping_key = "group",
#   motif_length = 2,
#   sequence_type = "nt",
#   output_type = "fraction"
# )
# print(result) # Returns a matrix of k-mer fractions per group
