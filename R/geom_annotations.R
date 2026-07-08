## -----------------------------------------------------------------
## msaVariant: dedicated annotation geoms (Zenodo-backed)
## -----------------------------------------------------------------
##
## Each geom accepts either `data` (a user-supplied data.frame) or
## `gene` (an HGNC symbol; the package fetches from the Zenodo
## deposit and uses the appropriate slice).
##
## For per-substitution sources (AlphaMissense, REVEL, CADD), the
## display defaults to a per-residue summary (mean across the 19
## possible alts at each position). The full per-substitution table
## is still available to users via the get_*() functions.

# Helper: turn a per-substitution data.frame into a per-residue
# tidy table for the heatmap-strip display. The summarisation
# choice is documented per geom.
.summarise_per_residue <- function(d, value_col, summary = "mean") {
  if (nrow(d) == 0L) return(d)
  agg_fn <- switch(summary,
                   "mean" = function(x) mean(x, na.rm = TRUE),
                   "max"  = function(x) max(x, na.rm = TRUE),
                   "median" = function(x) stats::median(x, na.rm = TRUE),
                   stop("Unknown summary method: ", summary))
  by_pos <- split(d[[value_col]], d$pos)
  out <- data.frame(
    pos = as.integer(names(by_pos)),
    stringsAsFactors = FALSE
  )
  out[[value_col]] <- vapply(by_pos, agg_fn, numeric(1))
  out
}

# Resolver: turns (`data`, `gene`) into a data.frame, fetching if
# the user passed `gene` and erroring on the wrong combinations.
.resolve_annotation <- function(data, gene, slot, geom_name) {
  if (!is.null(data) && !is.null(gene)) {
    rlang::abort(sprintf(
      "%s: pass either `data` or `gene`, not both.", geom_name))
  }
  if (is.null(data) && is.null(gene)) {
    rlang::abort(sprintf(
      "%s: must supply either `data` (a tidy data.frame) or `gene` (HGNC symbol).",
      geom_name))
  }
  if (!is.null(gene)) {
    bundle <- fetch_gene_data(gene)
    if (is.null(bundle)) {
      rlang::warn(sprintf(
        "%s: no data fetched for gene '%s'; layer will be empty.",
        geom_name, gene))
      return(NULL)
    }
    return(bundle[[slot]])
  }
  data
}

# ----- gnomAD ------------------------------------------------------

#' Overlay gnomAD allele frequencies onto an MSA
#'
#' Renders a per-residue heatmap strip of gnomAD allele frequencies.
#' For residues with multiple variants, the maximum joint allele
#' frequency across variants is plotted.
#'
#' @param data Optional pre-built data.frame; see
#'   `DATA_FORMAT_SPEC.md` for the schema. If `NULL`, supply `gene`.
#' @param gene HGNC symbol; triggers fetch from Zenodo deposit.
#' @param msa,ref_name Reference MSA and sequence name.
#' @param y_offset,track_height Track geometry.
#' @return A list of ggplot2 layers (or `NULL`).
#' @export
geom_gnomad <- function(data = NULL, gene = NULL,
                        msa, ref_name = NULL,
                        y_offset = -2, track_height = 1.2) {
  d <- .resolve_annotation(data, gene, "gnomad", "geom_gnomad")
  if (is.null(d) || nrow(d) == 0L) return(NULL)
  # Per-residue summary: max AF at each position
  if ("af_joint" %in% names(d)) {
    summary <- .summarise_per_residue(d, "af_joint", summary = "max")
    names(summary)[2] <- "max_af"
  } else {
    summary <- d
  }
  geom_track(summary, msa = msa, ref_name = ref_name,
             value = "max_af", type = "continuous",
             name = "gnomAD AF",
             palette = c("#FFFFFF", "#DADAEB", "#9E9AC8", "#6A51A3", "#3F007D"),
             y_offset = y_offset, track_height = track_height)
}

# ----- ClinVar -----------------------------------------------------

#' Overlay ClinVar pathogenicity calls onto an MSA
#'
#' Renders a discrete tick row showing ClinVar calls at each
#' affected residue. Where multiple ClinVar entries exist at one
#' residue, the most pathogenic call is shown.
#'
#' @inheritParams geom_gnomad
#' @param significance Character vector of significance levels to
#'   show (default: all five ACMG tiers).
#' @return A list of ggplot2 layers (or `NULL`).
#' @export
geom_clinvar <- function(data = NULL, gene = NULL,
                         msa, ref_name = NULL,
                         significance = c("Pathogenic",
                                          "Likely_pathogenic",
                                          "VUS",
                                          "Likely_benign",
                                          "Benign"),
                         y_offset = -3.5, track_height = 1.0) {
  d <- .resolve_annotation(data, gene, "clinvar", "geom_clinvar")
  if (is.null(d) || nrow(d) == 0L) return(NULL)
  d <- d[as.character(d$significance) %in% significance, , drop = FALSE]
  if (nrow(d) == 0L) {
    rlang::warn("geom_clinvar: no entries matched the requested significance levels.")
    return(NULL)
  }
  # If multiple entries at same residue, keep the most pathogenic
  sig_order <- c("Pathogenic" = 1, "Likely_pathogenic" = 2,
                 "VUS" = 3, "Likely_benign" = 4, "Benign" = 5)
  d$.priority <- sig_order[as.character(d$significance)]
  d <- d[order(d$.priority), ]
  d <- d[!duplicated(d$pos), ]
  d$.priority <- NULL

  palette <- c("Pathogenic"        = "#D7301F",
               "Likely_pathogenic" = "#FC8D59",
               "VUS"               = "#BDBDBD",
               "Likely_benign"     = "#91BFDB",
               "Benign"            = "#4575B4")
  geom_track(d, msa = msa, ref_name = ref_name,
             value = "significance", type = "discrete",
             name = "ClinVar", palette = palette,
             y_offset = y_offset, track_height = track_height)
}

# ----- AlphaMissense -----------------------------------------------

#' Overlay AlphaMissense per-residue mean scores onto an MSA
#'
#' Aggregates the per-substitution AlphaMissense table to a
#' per-residue mean for the heatmap strip. To inspect specific
#' substitutions, fetch the full table with `get_alphamissense()`.
#'
#' Note: AlphaMissense is licensed CC-BY-NC-SA 4.0.
#'
#' @inheritParams geom_gnomad
#' @param summary How to aggregate the 19 alt scores at each
#'   position. One of `"mean"` (default), `"max"`, `"median"`.
#' @return A list of ggplot2 layers (or `NULL`).
#' @export
geom_alphamissense <- function(data = NULL, gene = NULL,
                               msa, ref_name = NULL,
                               summary = c("mean", "max", "median"),
                               y_offset = -5, track_height = 1.2) {
  summary <- match.arg(summary)
  d <- .resolve_annotation(data, gene, "alphamissense", "geom_alphamissense")
  if (is.null(d) || nrow(d) == 0L) return(NULL)
  d2 <- .summarise_per_residue(d, "am_score", summary = summary)
  names(d2)[2] <- "score"
  geom_track(d2, msa = msa, ref_name = ref_name,
             value = "score", type = "continuous",
             name = sprintf("AlphaMissense (%s)", summary),
             value_range = c(0, 1),
             palette = c("#2c7bb6","#abd9e9","#ffffbf","#fdae61","#d7191c"),
             y_offset = y_offset, track_height = track_height)
}

# ----- REVEL -------------------------------------------------------

#' Overlay REVEL per-residue scores onto an MSA
#'
#' @inheritParams geom_alphamissense
#' @return A list of ggplot2 layers (or `NULL`).
#' @export
geom_revel <- function(data = NULL, gene = NULL,
                       msa, ref_name = NULL,
                       summary = c("mean", "max", "median"),
                       y_offset = -6.5, track_height = 1.2) {
  summary <- match.arg(summary)
  d <- .resolve_annotation(data, gene, "revel", "geom_revel")
  if (is.null(d) || nrow(d) == 0L) return(NULL)
  d2 <- .summarise_per_residue(d, "revel_score", summary = summary)
  names(d2)[2] <- "score"
  geom_track(d2, msa = msa, ref_name = ref_name,
             value = "score", type = "continuous",
             name = sprintf("REVEL (%s)", summary),
             value_range = c(0, 1),
             palette = c("#2c7bb6","#abd9e9","#ffffbf","#fdae61","#d7191c"),
             y_offset = y_offset, track_height = track_height)
}

# ----- CADD --------------------------------------------------------

#' Overlay CADD per-residue PHRED scores onto an MSA
#'
#' CADD is licensed for non-commercial use only.
#'
#' @inheritParams geom_alphamissense
#' @param value_range Length-2 PHRED range for the colour scale.
#'   Default `c(0, 40)` covers the practically informative range.
#' @return A list of ggplot2 layers (or `NULL`).
#' @export
geom_cadd <- function(data = NULL, gene = NULL,
                      msa, ref_name = NULL,
                      summary = c("median", "mean", "max"),
                      value_range = c(0, 40),
                      y_offset = -8, track_height = 1.2) {
  summary <- match.arg(summary)
  d <- .resolve_annotation(data, gene, "cadd", "geom_cadd")
  if (is.null(d) || nrow(d) == 0L) return(NULL)
  d2 <- .summarise_per_residue(d, "cadd_phred", summary = summary)
  names(d2)[2] <- "score"
  geom_track(d2, msa = msa, ref_name = ref_name,
             value = "score", type = "continuous",
             name = sprintf("CADD PHRED (%s)", summary),
             value_range = value_range,
             palette = c("#f7fcb9","#addd8e","#41ab5d","#238443","#005a32"),
             y_offset = y_offset, track_height = track_height)
}
