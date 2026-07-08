## -----------------------------------------------------------------
## msaVariant: coordinate mapping
## -----------------------------------------------------------------
##
## A clinical variant is annotated against a reference protein (e.g.
## "PATL1 p.K518fs" against ENSP00000300146). The MSA, however, has
## its own column numbering and contains gaps. Every annotation
## function in this package therefore needs to translate
## (reference residue position) -> (MSA column).
##
## This module provides that translation and helpers around it.

#' Build a residue -> MSA column lookup
#'
#' For one row of a multiple sequence alignment (the "reference"
#' sequence), build a map from ungapped residue position
#' (1, 2, 3, ...) to alignment column.
#'
#' @param msa An object understood by ggmsa: an AAMultipleAlignment,
#'   DNAMultipleAlignment, AAStringSet, DNAStringSet, character vector
#'   of equal-length sequences, or a path to a FASTA file.
#' @param ref_name Name of the reference sequence in the MSA. The
#'   default `NULL` uses the first sequence.
#' @return A data.frame with columns `residue_pos`, `aa`, `msa_col`.
#'   Gap columns in the reference are skipped.
#' @export
build_msa_coord_map <- function(msa, ref_name = NULL) {
  seqs <- .coerce_to_named_char(msa)

  if (is.null(ref_name)) {
    ref_name <- names(seqs)[1L]
  }
  if (!ref_name %in% names(seqs)) {
    rlang::abort(sprintf(
      "Reference '%s' not found in MSA. Available: %s",
      ref_name, paste(names(seqs), collapse = ", ")
    ))
  }
  ref_seq <- seqs[[ref_name]]
  chars   <- strsplit(ref_seq, "", fixed = TRUE)[[1L]]

  is_gap  <- chars %in% c("-", ".", "*")
  ungapped_pos <- cumsum(!is_gap)
  ungapped_pos[is_gap] <- NA_integer_

  data.frame(
    residue_pos = ungapped_pos,
    aa          = chars,
    msa_col     = seq_along(chars),
    stringsAsFactors = FALSE
  )
}

#' Map a vector of variant residue positions to MSA columns
#'
#' Thin wrapper around `build_msa_coord_map()` that vectorises the
#' lookup and warns when a position is out of range.
#'
#' @param positions Integer vector of reference protein positions.
#' @param msa,ref_name See `build_msa_coord_map()`.
#' @return Integer vector of MSA column indices (NA where mapping
#'   fails).
#' @export
map_variant_to_msa <- function(positions, msa, ref_name = NULL) {
  m <- build_msa_coord_map(msa, ref_name)
  # invert: residue_pos -> msa_col, drop NA residue rows (gap)
  m <- m[!is.na(m$residue_pos), , drop = FALSE]
  lookup <- stats::setNames(m$msa_col, as.character(m$residue_pos))
  out <- unname(lookup[as.character(positions)])

  n_miss <- sum(is.na(out) & !is.na(positions))
  if (n_miss > 0L) {
    rlang::warn(sprintf(
      "%d position(s) could not be mapped to the reference MSA (out of range or gap).",
      n_miss
    ))
  }
  out
}

# Internal helper: take whatever `msa` representation is supplied
# and return a named character vector of equal-length sequences.
.coerce_to_named_char <- function(msa) {
  if (is.character(msa) && length(msa) == 1L && file.exists(msa)) {
    # FASTA path
    aa <- Biostrings::readAAStringSet(msa)
    seqs <- as.character(aa)
    if (is.null(names(seqs)) || any(names(seqs) == "")) {
      names(seqs) <- paste0("seq", seq_along(seqs))
    }
    return(seqs)
  }
  if (inherits(msa, c("AAMultipleAlignment", "DNAMultipleAlignment"))) {
    return(as.character(msa@unmasked))
  }
  if (inherits(msa, c("AAStringSet", "DNAStringSet", "BStringSet"))) {
    return(as.character(msa))
  }
  if (is.character(msa)) {
    if (is.null(names(msa)) || any(names(msa) == "")) {
      names(msa) <- paste0("seq", seq_along(msa))
    }
    return(msa)
  }
  rlang::abort(sprintf(
    "Don't know how to extract sequences from object of class '%s'",
    paste(class(msa), collapse = "/")
  ))
}
