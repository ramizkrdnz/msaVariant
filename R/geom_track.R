## -----------------------------------------------------------------
## msaVariant: geom_track()
## -----------------------------------------------------------------
##
## A single generic layer that plots any per-residue annotation
## the user brings themselves -- continuous (e.g. AlphaMissense
## scores, gnomAD AF) or categorical (e.g. ClinVar significance).
##
## Replaces the geom_gnomad / geom_clinvar / geom_alphamissense
## / geom_revel / geom_cadd layers that v0.1 deliberately does NOT
## ship, because we cannot ship reliable per-gene snapshots of
## those databases inside the package.
##
## Example usage:
##   my_gnomad <- read.csv("my_gnomad_query.csv")  # cols: pos, af
##   geom_track(my_gnomad, msa = fa, value = "af",
##              name = "gnomAD AF", type = "continuous")

#' Overlay an arbitrary per-residue annotation track onto an MSA
#'
#' Renders a strip below (or above) the MSA showing any per-residue
#' scalar or categorical value the user supplies. Handles both
#' continuous data (rendered as a colour-graded heatmap strip) and
#' discrete data (rendered as a coloured tick row).
#'
#' This is the generic, bring-your-own-data layer for annotation
#' tracks. Use it for gnomAD allele frequencies, ClinVar significance,
#' AlphaMissense / REVEL / CADD scores, post-translational-modification
#' sites, ChIP-seq signal, or anything else you can map to a residue.
#'
#' @param data A data.frame with at least:
#'   * `pos`   — integer residue position (1-based in the reference)
#'   * a value column whose name you pass via `value`
#' @param msa,ref_name MSA and reference sequence name (same as in
#'   the parent `ggmsa()` call).
#' @param value Name of the column in `data` holding the
#'   per-residue value to plot.
#' @param type Either `"continuous"` (heatmap strip, default) or
#'   `"discrete"` (categorical tick row).
#' @param name Label printed to the right of the track.
#' @param palette For `type = "continuous"`, a vector of colours
#'   passed to `colorRampPalette()`. For `type = "discrete"`, a
#'   named vector mapping levels of `value` to colours.
#' @param value_range For continuous tracks, a length-2 numeric vector
#'   giving `c(min, max)` for the colour scale. Defaults to the
#'   observed range.
#' @param y_offset,track_height Vertical placement and height in
#'   MSA-row units.
#' @return A list of ggplot2 layers.
#' @export
#' @examples
#' \dontrun{
#' # Continuous: gnomAD allele frequencies from your own query
#' my_af <- data.frame(pos = 510:530, af = runif(21, 0, 1e-3))
#' geom_track(my_af, msa = fa, value = "af",
#'            name = "gnomAD AF", type = "continuous")
#'
#' # Discrete: ClinVar calls from your own export
#' my_cv <- data.frame(pos = c(498, 518),
#'                     sig = c("Benign", "Pathogenic"))
#' geom_track(my_cv, msa = fa, value = "sig", type = "discrete",
#'            name = "ClinVar",
#'            palette = c(Benign = "#4575B4", Pathogenic = "#D7301F"))
#' }
geom_track <- function(data,
                       msa,
                       ref_name     = NULL,
                       value,
                       type         = c("continuous", "discrete"),
                       name         = "",
                       palette      = NULL,
                       value_range  = NULL,
                       y_offset     = -2,
                       track_height = 1.2) {
  type <- match.arg(type)
  if (!is.data.frame(data)) {
    rlang::abort("`data` must be a data.frame.")
  }
  if (!"pos" %in% names(data)) {
    rlang::abort("`data` must contain a `pos` column.")
  }
  if (missing(value) || !value %in% names(data)) {
    rlang::abort(sprintf("Column '%s' not found in `data`.", value))
  }

  # Map residue position to MSA column
  data$msa_col <- map_variant_to_msa(data$pos, msa, ref_name)
  data <- data[!is.na(data$msa_col), , drop = FALSE]
  if (nrow(data) == 0L) {
    rlang::warn("No track positions mapped to MSA columns.")
    return(NULL)
  }

  y0 <- y_offset
  y1 <- y_offset + track_height

  vals <- data[[value]]

  if (type == "continuous") {
    if (!is.numeric(vals)) {
      rlang::abort("Continuous track requires a numeric value column.")
    }
    if (is.null(palette)) {
      palette <- c("#2c7bb6","#abd9e9","#ffffbf","#fdae61","#d7191c")
    }
    if (is.null(value_range)) value_range <- range(vals, na.rm = TRUE)
    if (diff(value_range) == 0) value_range <- value_range + c(-0.5, 0.5)
    clipped <- pmin(pmax(vals, value_range[1]), value_range[2])
    norm    <- (clipped - value_range[1]) /
               (value_range[2] - value_range[1])
    fills   <- grDevices::colorRampPalette(palette)(101)[
                 pmax(1, pmin(101, round(norm * 100) + 1))]
  } else {
    # discrete
    if (is.null(palette)) {
      lv <- unique(as.character(vals))
      palette <- setNames(
        scales::hue_pal()(length(lv)),
        lv
      )
    }
    fills <- unname(palette[as.character(vals)])
    fills[is.na(fills)] <- "#888888"
  }

  tiles <- ggplot2::geom_rect(
    data = data.frame(
      xmin = data$msa_col - 0.5,
      xmax = data$msa_col + 0.5,
      ymin = y0, ymax = y1,
      stringsAsFactors = FALSE
    ),
    mapping = ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                           ymin = .data$ymin, ymax = .data$ymax),
    inherit.aes = FALSE,
    fill        = fills
  )

  label <- ggplot2::annotate(
    "text",
    x = max(data$msa_col) + 1.5,
    y = (y0 + y1) / 2,
    label = name,
    hjust = 0, size = 3.2, fontface = "italic"
  )

  list(tiles, label)
}
