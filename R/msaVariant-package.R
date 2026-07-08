#' msaVariant: Clinical-genetics MSA visualisation with variant overlay
#'
#' Extends `ggmsa` with annotation layers tailored to clinical and
#' rare-disease genetics workflows.
#'
#' What the package provides (v0.1):
#'
#' * `geom_variant()`        — user variants (tidy data frame)
#' * `geom_domain()`         — protein domain track (user-supplied)
#' * `geom_track()`          — generic per-residue annotation track
#'                             for any user-supplied data (gnomAD AF,
#'                             ClinVar significance, AlphaMissense
#'                             scores, PTM sites, etc.)
#' * `conservation_score()`  — per-column Shannon / Jensen-Shannon
#'                             computed from the MSA itself
#'
#' Plus utility functions:
#'
#' * `build_msa_coord_map()` — residue-to-column lookup
#' * `map_variant_to_msa()`  — vectorised position mapping
#' * `scale_fill_pathogenicity()` — colourblind-safe ACMG palette
#'
#' What the package does NOT bundle:
#'
#' Population-genetics and pathogenicity databases (gnomAD, ClinVar,
#' AlphaMissense, REVEL, CADD) change frequently, are licensed in
#' ways that constrain redistribution, and are too large to ship as
#' static files. Users supply their own annotation data via
#' `geom_track()`. A future v0.2 will add optional online fetchers
#' that pull these databases on demand and cache results locally.
#'
#' @name msaVariant-package
#' @aliases msaVariant
#' @importFrom stats aggregate median setNames
"_PACKAGE"

## Silence R CMD check NOTEs about non-standard evaluation: these names
## are ggplot2 `aes()` mappings and the rlang `.data` pronoun, referenced
## as bare symbols rather than defined as top-level bindings.
utils::globalVariables(c(
  ".data", "msa_col", "seq_name", "residue", "x", "y",
  "code", "txt", "name", "msa_start", "msa_end"
))
