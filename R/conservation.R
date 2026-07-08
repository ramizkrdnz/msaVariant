## -----------------------------------------------------------------
## msaVariant: conservation scoring
## -----------------------------------------------------------------
##
## Per-column conservation is essential context for variant
## interpretation: a variant at a Shannon-conserved column carries
## stronger functional evidence than one at a variable position.
##
## Two metrics implemented here:
##   * Shannon entropy (the standard, easy to interpret)
##   * Jensen-Shannon divergence vs a background distribution
##     (slightly better for proteins; corrects for amino-acid abundance)

#' Per-column conservation score for an MSA
#'
#' @param msa MSA in any form accepted by `build_msa_coord_map()`.
#' @param method One of `"shannon"` (default) or `"js"`
#'   (Jensen-Shannon divergence against a background).
#' @param background Optional named vector of background residue
#'   frequencies. If `NULL` and `method = "js"`, a uniform background
#'   is used.
#' @param ignore_gaps If `TRUE` (default), gap characters are
#'   excluded from the per-column frequency calculation.
#' @return A data.frame with columns `msa_col` and `score`.
#'   For Shannon entropy the score is in nats, rescaled so 1 = fully
#'   conserved and 0 = uniform across the observed alphabet. For
#'   Jensen-Shannon, score is in [0, 1] with 1 = fully conserved.
#' @export
conservation_score <- function(msa,
                               method      = c("shannon", "js"),
                               background  = NULL,
                               ignore_gaps = TRUE) {
  method <- match.arg(method)
  seqs   <- .coerce_to_named_char(msa)

  mat <- do.call(rbind, strsplit(seqs, "", fixed = TRUE))
  if (length(unique(nchar(seqs))) != 1L) {
    rlang::abort("Sequences in MSA are not equal length.")
  }

  ncol_msa <- ncol(mat)
  score <- numeric(ncol_msa)

  for (j in seq_len(ncol_msa)) {
    col <- mat[, j]
    if (ignore_gaps) col <- col[!col %in% c("-", ".", "*")]
    if (length(col) == 0L) { score[j] <- NA_real_; next }

    freq <- table(col) / length(col)
    score[j] <- switch(
      method,
      shannon = .shannon_conservation(freq),
      js      = .js_conservation(freq, background)
    )
  }

  data.frame(msa_col = seq_len(ncol_msa),
             score   = score,
             stringsAsFactors = FALSE)
}

# Shannon entropy converted to a conservation score in [0, 1].
# H_max is computed against the full canonical alphabet size (20 for
# proteins, 4 for nucleotides) -- NOT the observed alphabet -- so
# that columns dominated by one residue with one or two outliers
# still score as highly conserved. Detection is heuristic from the
# observed residues; user can override with `alphabet_size`.
.shannon_conservation <- function(freq, alphabet_size = NULL) {
  if (length(freq) <= 1L) return(1.0)
  if (is.null(alphabet_size)) {
    aa_set <- c("A","C","D","E","F","G","H","I","K","L",
                "M","N","P","Q","R","S","T","V","W","Y")
    is_protein <- all(names(freq) %in% c(aa_set, "X", "B", "Z", "U", "O"))
    alphabet_size <- if (is_protein) 20L else 4L
  }
  H     <- -sum(freq * log(freq))
  H_max <- log(alphabet_size)
  max(0, 1.0 - H / H_max)
}

# Jensen-Shannon divergence against a background distribution,
# normalised to [0, 1].
.js_conservation <- function(freq, background = NULL) {
  alphabet <- names(freq)
  if (is.null(background)) {
    background <- rep(1 / length(alphabet), length(alphabet))
    names(background) <- alphabet
  } else {
    background <- background[alphabet]
    background[is.na(background)] <- 1e-9
    background <- background / sum(background)
  }
  p <- as.numeric(freq)
  q <- as.numeric(background)
  m <- 0.5 * (p + q)
  kl <- function(a, b) sum(ifelse(a > 0, a * log2(a / b), 0))
  js <- 0.5 * kl(p, m) + 0.5 * kl(q, m)
  # JS in [0, 1] with log2 base when comparing two distributions
  js
}
