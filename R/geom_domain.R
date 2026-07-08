## -----------------------------------------------------------------
## msaVariant: geom_domain()
## -----------------------------------------------------------------

#' Overlay protein domain annotation on an MSA
#'
#' Draws a horizontal track showing protein domain ranges with
#' labels, mapped to MSA columns via the reference sequence.
#'
#' @param domains data.frame with columns `start`, `end`, `name`
#'   (residue positions, 1-based, inclusive). Optional column
#'   `source` (e.g. "InterPro", "Pfam", "SMART") is preserved.
#'   If `NULL`, supply `gene` instead.
#' @param gene HGNC gene symbol; if supplied, the package will
#'   fetch domain annotations from the Zenodo deposit.
#' @param msa,ref_name Reference MSA and sequence.
#' @param y_offset,track_height Track geometry.
#' @param fill Default fill colour for domain boxes if no
#'   `domains$colour` column is provided.
#' @return A list of ggplot2 layers.
#' @export
geom_domain <- function(domains = NULL,
                        gene    = NULL,
                        msa,
                        ref_name     = NULL,
                        y_offset     = NULL,
                        track_height = 1.0,
                        fill         = "#999999") {
  if (is.null(domains) && is.null(gene)) {
    rlang::abort("geom_domain: supply either `domains` or `gene`.")
  }
  if (!is.null(domains) && !is.null(gene)) {
    rlang::abort("geom_domain: supply `domains` or `gene`, not both.")
  }
  if (!is.null(gene)) {
    bundle <- fetch_gene_data(gene)
    if (is.null(bundle)) {
      rlang::warn(sprintf("geom_domain: no data bundle for '%s'.", gene))
      return(NULL)
    }
    domains <- bundle$domains
    if (is.null(domains) || nrow(domains) == 0L) {
      rlang::warn(sprintf("geom_domain: no domains in bundle for '%s'.", gene))
      return(NULL)
    }
  }

  required <- c("start", "end", "name")
  missing  <- setdiff(required, names(domains))
  if (length(missing)) {
    rlang::abort(sprintf("domains missing column(s): %s",
                         paste(missing, collapse = ", ")))
  }

  d <- domains
  d$msa_start <- map_variant_to_msa(d$start, msa, ref_name)
  d$msa_end   <- map_variant_to_msa(d$end,   msa, ref_name)
  d <- d[!is.na(d$msa_start) & !is.na(d$msa_end), , drop = FALSE]
  if (nrow(d) == 0L) {
    rlang::warn("No domains mapped to MSA columns.")
    return(NULL)
  }

  seqs <- .coerce_to_named_char(msa)
  n_seq <- length(seqs)
  y0 <- if (is.null(y_offset)) n_seq + 1.7 else y_offset
  y1 <- y0 + track_height

  fills <- if ("colour" %in% names(d)) d$colour else rep(fill, nrow(d))

  boxes <- ggplot2::geom_rect(
    data = data.frame(
      xmin = d$msa_start - 0.5,
      xmax = d$msa_end   + 0.5,
      ymin = y0, ymax = y1,
      fill = fills,
      stringsAsFactors = FALSE
    ),
    mapping = ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                           ymin = .data$ymin, ymax = .data$ymax),
    inherit.aes = FALSE,
    fill        = fills,
    colour      = "black", linewidth = 0.3
  )

  labels <- ggplot2::annotate(
    "text",
    x = (d$msa_start + d$msa_end) / 2,
    y = (y0 + y1) / 2,
    label = d$name,
    size  = 3.3, fontface = "bold", colour = "white"
  )

  list(boxes, labels)
}
