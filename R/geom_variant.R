## -----------------------------------------------------------------
## msaVariant: geom_variant()
## -----------------------------------------------------------------
##
## The flagship layer. Given a data frame of variants (or a VCF + a
## transcript), draws markers at the corresponding MSA columns.
##
## Design notes:
## * We render the marker as a vertical band that spans the entire
##   alignment (so it's visible whether the user is plotting 5 or
##   500 sequences) plus a labeled tick at the top.
## * Pathogenicity classification controls colour. Default colour
##   scale is colorblind-safe Okabe-Ito-derived.
## * For frameshifts/indels we draw a span (start..end of affected
##   residues) rather than a point.

#' Overlay variants onto a multiple sequence alignment
#'
#' Adds vertical bands marking the alignment columns affected by one
#' or more clinical / population variants. Works with a tidy data
#' frame; for VCF input use `read_vcf_variants()` first.
#'
#' @param variants A data.frame with at minimum a column `pos`
#'   (integer, 1-based protein residue position in the reference).
#'   Optional columns:
#'   * `pos_end`   : integer, for indels/frameshifts; defaults to
#'                   `pos` (single-residue variant)
#'   * `label`     : character, shown above the band
#'   * `consequence` : factor, used to colour the band
#'                     (e.g. "missense", "frameshift", "nonsense",
#'                     "synonymous", "splice")
#'   * `pathogenicity` : factor, alternative colour key (e.g.
#'                     "Pathogenic", "Likely_pathogenic", "VUS",
#'                     "Likely_benign", "Benign")
#' @param msa,ref_name The MSA used in the parent `ggmsa()` call,
#'   and the reference sequence within it. The geom maps variant
#'   positions to alignment columns via `map_variant_to_msa()`.
#' @param colour_by Either `"consequence"` (default) or
#'   `"pathogenicity"`, controlling which column drives the fill.
#' @param alpha Band transparency in `[0, 1]`. Default `0.5`.
#' @param show_label Logical; whether to draw labels above the band.
#' @param label_y Numeric vertical position of the variant label.
#'   Defaults to `n_seq + 3` (just above the MSA, clear of the
#'   default `geom_domain()` track at `n_seq + 1.7`).
#' @param ... Passed to the underlying `geom_rect()` layer.
#' @return A ggplot2 layer (or list of layers) that can be `+`-ed
#'   onto a `ggmsa()` plot.
#' @export
#' @examples
#' \dontrun{
#' library(ggmsa)
#' fa <- system.file("extdata", "patl1_orthologs.fasta",
#'                   package = "msaVariant")
#' v  <- data.frame(pos = 518, pos_end = 577,
#'                  label = "p.K518fs",
#'                  consequence = "frameshift")
#' ggmsa(fa, start = 500, end = 580) +
#'   geom_variant(v, msa = fa, ref_name = "PATL1_HUMAN")
#' }
geom_variant <- function(variants,
                         msa,
                         ref_name   = NULL,
                         colour_by  = c("consequence", "pathogenicity"),
                         alpha      = 0.5,
                         show_label = TRUE,
                         label_y    = NULL,
                         ...) {
  colour_by <- match.arg(colour_by)
  if (!is.data.frame(variants)) {
    rlang::abort("`variants` must be a data.frame.")
  }
  if (!"pos" %in% names(variants)) {
    rlang::abort("`variants` must contain a `pos` column (integer residue position).")
  }
  if (!"pos_end" %in% names(variants)) variants$pos_end <- variants$pos
  if (!"label"   %in% names(variants)) variants$label   <- ""

  variants$msa_xmin <- map_variant_to_msa(variants$pos,     msa, ref_name) - 0.5
  variants$msa_xmax <- map_variant_to_msa(variants$pos_end, msa, ref_name) + 0.5

  variants <- variants[!is.na(variants$msa_xmin) & !is.na(variants$msa_xmax), ]
  if (nrow(variants) == 0L) {
    rlang::warn("No variants mapped to MSA columns; geom_variant() will be empty.")
    return(NULL)
  }

  seqs  <- .coerce_to_named_char(msa)
  n_seq <- length(seqs)

  # Resolve fill colour via the shared palette, hard-coding it onto
  # the layer so this geom never collides with another geom's fill
  # scale. Users wanting a different palette can override via
  # `variants$colour`.
  palette_lookup <- c(
    "Pathogenic" = "#D7301F", "Likely_pathogenic" = "#FC8D59",
    "VUS" = "#BDBDBD", "Likely_benign" = "#91BFDB", "Benign" = "#4575B4",
    "frameshift" = "#7F3B08", "nonsense" = "#B35806",
    "missense"   = "#F1A340", "synonymous" = "#998EC3",
    "splice"     = "#542788"
  )
  key <- if (colour_by %in% names(variants)) as.character(variants[[colour_by]]) else NA_character_
  fills_resolved <- if ("colour" %in% names(variants)) {
    variants$colour
  } else {
    unname(palette_lookup[key])
  }
  fills_resolved[is.na(fills_resolved)] <- "#888888"

  # The band must span the full plot panel regardless of whether
  # the y-axis is continuous (geom_tile-based) or discrete (ggmsa).
  # Using ymin=-Inf, ymax=Inf is the ggplot2 idiom for "fill the
  # panel height" -- it works on both discrete and continuous y
  # scales.
  rect_layer <- ggplot2::geom_rect(
    data = data.frame(
      xmin = variants$msa_xmin,
      xmax = variants$msa_xmax,
      stringsAsFactors = FALSE
    ),
    mapping = ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                           ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill        = fills_resolved,
    alpha       = alpha,
    ...
  )

  layers <- list(rect_layer)

  if (show_label && any(nzchar(variants$label))) {
    # Place the label above the MSA panel. We use a numeric
    # y-offset above n_seq (same convention as geom_domain),
    # which `ggmsa` correctly renders above the alignment.
    # The default n_seq + 3 keeps the label clear of the
    # geom_domain() track at n_seq + 1.7.
    y_label <- if (is.null(label_y)) n_seq + 3 else label_y
    label_layer <- ggplot2::annotate(
      "text",
      x      = (variants$msa_xmin + variants$msa_xmax) / 2,
      y      = y_label,
      label  = variants$label,
      size   = 4.0,
      fontface = "bold",
      vjust  = 0
    )
    layers <- c(layers, list(label_layer))
  }

  layers
}
