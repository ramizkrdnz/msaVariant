## -----------------------------------------------------------------
## msaVariant: shared scales
## -----------------------------------------------------------------

#' Default colour scale for ACMG pathogenicity tiers
#'
#' Provides a colourblind-safe 5-tier palette ordered from
#' pathogenic (red) to benign (blue), with grey for VUS.
#'
#' @param ... Passed to `ggplot2::scale_fill_manual()`.
#' @return A `ScaleDiscrete` object.
#' @examples
#' library(ggplot2)
#' d <- data.frame(x = c("Pathogenic", "VUS", "Benign"), y = 1)
#' ggplot(d, aes(x, y, fill = x)) +
#'   geom_col() +
#'   scale_fill_pathogenicity()
#' @export
scale_fill_pathogenicity <- function(...) {
  ggplot2::scale_fill_manual(
    name   = "ClinVar significance",
    values = c(
      "Pathogenic"          = "#D7301F",
      "Likely_pathogenic"   = "#FC8D59",
      "VUS"                 = "#BDBDBD",
      "Likely_benign"       = "#91BFDB",
      "Benign"              = "#4575B4",
      "constrained"         = "#7B3294",
      "observed"            = "#C2A5CF",
      "frameshift"          = "#7F3B08",
      "nonsense"            = "#B35806",
      "missense"            = "#F1A340",
      "synonymous"          = "#998EC3",
      "splice"              = "#542788"
    ),
    drop = FALSE,
    ...
  )
}
