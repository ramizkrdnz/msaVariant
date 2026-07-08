## ===================================================================
## msaVariant: compute_acmg_codes()
## ===================================================================
##
## Pure, testable evaluation of a small set of ACMG/AMP evidence codes
## for a single amino-acid substitution, computed ENTIRELY from the
## tables already present in a gene bundle (meta/domains/clinvar/
## gnomad/alphamissense/revel/cadd). No new data source is introduced.
##
## The rules implemented here are the exact operational definitions
## requested for the annotation strip; they are intentionally simple
## and transparent, not a full ACMG classifier.

#' Compute triggered ACMG evidence codes for one variant
#'
#' Evaluates a fixed subset of ACMG/AMP evidence codes for a single
#' missense substitution using only the annotation tables carried in a
#' gene bundle. Returns just the codes that are triggered.
#'
#' Rules (all computed from the bundle; no external lookups):
#' \describe{
#'   \item{PS1}{ClinVar contains the exact same \code{aa_change} with
#'     significance \code{"Pathogenic"}.}
#'   \item{PM1}{\code{variant_pos} falls within the start–end range of
#'     at least one \code{domains} row.}
#'   \item{PM2}{This \code{aa_change} is absent from gnomAD, or its
#'     \code{af_joint} is below 0.0001.}
#'   \item{PM5}{ClinVar has a record at the same position with a
#'     \emph{different} alt residue and significance
#'     \code{"Pathogenic"} or \code{"Likely_pathogenic"}.}
#'   \item{PP3}{At this substitution, count how many of the three
#'     computational predictors pass their threshold — AlphaMissense
#'     \code{am_class == "likely_pathogenic"}, REVEL
#'     \code{revel_score > 0.7}, CADD \code{cadd_phred >= 25}. Triggered
#'     when the passing count is at least \code{pp3_min_predictors}.}
#' }
#'
#' Codes whose evaluation requires the alternate allele (PS1, PM2, PM5,
#' PP3) are only evaluated when \code{aa_change} parses to a full
#' ref+pos+alt substitution (e.g. \code{"R175H"}). PM1 depends only on
#' position and is always evaluated. This keeps the function safe when
#' called with a position-only label (it simply won't over-fire PM2).
#'
#' @param bundle A gene bundle: a named list with elements
#'   \code{domains}, \code{clinvar}, \code{gnomad},
#'   \code{alphamissense}, \code{revel}, \code{cadd} (as produced by
#'   [import_local_bundle()] / [fetch_gene_data()]).
#' @param variant_pos Integer residue position (canonical UniProt
#'   numbering).
#' @param aa_change Character variant identifier, e.g. \code{"R175H"}
#'   or \code{"p.R175H"}. Ref/alt are parsed from this string.
#' @param pp3_min_predictors Integer 1–3 (default 2). PP3 fires when at
#'   least this many of the three PP3 predictors pass their threshold.
#'   \code{1} = any predictor, \code{2} = majority (default),
#'   \code{3} = all three must agree.
#'
#' @return A character vector of the triggered codes, in canonical
#'   evidence-strength order (\code{PS1}, \code{PM1}, \code{PM2},
#'   \code{PM5}, \code{PP3}); \code{character(0)} if none fire.
#'
#' @examples
#' \dontrun{
#' b <- fetch_gene_data("TP53")
#' compute_acmg_codes(b, 175, "R175H")
#' }
#'
#' @export
compute_acmg_codes <- function(bundle, variant_pos, aa_change,
                               pp3_min_predictors = 2L) {
  variant_pos <- as.integer(variant_pos)
  pp3_min_predictors <- as.integer(pp3_min_predictors)
  if (is.na(pp3_min_predictors) || pp3_min_predictors < 1L)
    pp3_min_predictors <- 1L
  if (pp3_min_predictors > 3L) pp3_min_predictors <- 3L

  ## --- Parse aa_change into ref / pos / alt --------------------------
  parsed <- .acmg_parse_aa_change(aa_change)
  ## If aa_change encodes a position, prefer the explicit variant_pos
  ## but fall back to the parsed one when variant_pos is missing.
  if (is.na(variant_pos) && !is.na(parsed$pos)) variant_pos <- parsed$pos
  has_alt <- !is.na(parsed$alt) && nzchar(parsed$alt)

  ## Canonical aa_change string used for exact matching. Prefer the
  ## caller's string (stripped of the "p." prefix); if only pos/alt are
  ## known, reconstruct it.
  aac <- sub("^p\\.", "", as.character(aa_change))

  triggered <- character(0)

  ## Small helpers -----------------------------------------------------
  nonempty <- function(df) !is.null(df) && is.data.frame(df) && nrow(df) > 0L
  has_cols <- function(df, cols) all(cols %in% names(df))
  chr <- function(x) as.character(x)

  ## --- PS1: exact same aa_change, ClinVar Pathogenic -----------------
  if (has_alt) {
    cv <- bundle$clinvar
    if (nonempty(cv) && has_cols(cv, c("aa_change", "significance"))) {
      hit <- chr(cv$aa_change) == aac & chr(cv$significance) == "Pathogenic"
      if (any(hit, na.rm = TRUE)) triggered <- c(triggered, "PS1")
    }
  }

  ## --- PM1: variant position within a domain range -------------------
  dom <- bundle$domains
  if (nonempty(dom) && has_cols(dom, c("start", "end")) && !is.na(variant_pos)) {
    in_dom <- dom$start <= variant_pos & dom$end >= variant_pos
    if (any(in_dom, na.rm = TRUE)) triggered <- c(triggered, "PM1")
  }

  ## --- PM2: absent from gnomAD OR af_joint < 1e-4 --------------------
  if (has_alt) {
    gn <- bundle$gnomad
    if (nonempty(gn) && has_cols(gn, c("aa_change", "af_joint"))) {
      rows <- gn[chr(gn$aa_change) == aac, , drop = FALSE]
      if (nrow(rows) == 0L) {
        triggered <- c(triggered, "PM2")                 # absent
      } else if (any(rows$af_joint < 1e-4, na.rm = TRUE)) {
        triggered <- c(triggered, "PM2")                 # rare
      }
    } else {
      ## No gnomAD table at all -> the variant is absent from gnomAD.
      triggered <- c(triggered, "PM2")
    }
  }

  ## --- PM5: same pos, different alt, ClinVar (Likely_)pathogenic -----
  if (has_alt) {
    cv <- bundle$clinvar
    if (nonempty(cv) && has_cols(cv, c("pos", "aa_alt", "significance")) &&
        !is.na(variant_pos)) {
      hit <- cv$pos == variant_pos &
             chr(cv$aa_alt) != parsed$alt &
             chr(cv$significance) %in% c("Pathogenic", "Likely_pathogenic")
      if (any(hit, na.rm = TRUE)) triggered <- c(triggered, "PM5")
    }
  }

  ## --- PP3: count deleterious predictors, trigger on >= threshold ---
  ## Predictors: AlphaMissense class == "likely_pathogenic",
  ## REVEL > 0.7, CADD phred >= 25. PP3 fires when at least
  ## pp3_min_predictors of them pass.
  if (has_alt) {
    am <- bundle$alphamissense
    rv <- bundle$revel
    cd <- bundle$cadd

    am_ok <- nonempty(am) && has_cols(am, c("aa_change", "am_class")) &&
      any(chr(am$aa_change) == aac &
          chr(am$am_class) == "likely_pathogenic", na.rm = TRUE)

    revel_ok <- nonempty(rv) && has_cols(rv, c("aa_change", "revel_score")) &&
      any(chr(rv$aa_change) == aac & rv$revel_score > 0.7, na.rm = TRUE)

    cadd_ok <- nonempty(cd) && has_cols(cd, c("aa_change", "cadd_phred")) &&
      any(chr(cd$aa_change) == aac & cd$cadd_phred >= 25, na.rm = TRUE)

    pp3_pass <- sum(c(am_ok, revel_ok, cadd_ok))
    if (pp3_pass >= pp3_min_predictors) triggered <- c(triggered, "PP3")
  }

  ## --- Return in canonical strength order ---------------------------
  order_ref <- c("PS1", "PM1", "PM2", "PM5", "PP3")
  triggered <- order_ref[order_ref %in% triggered]
  triggered
}

## Internal: parse "p.R175H" / "R175H" -> list(ref, pos, alt).
## Alt may be a single AA, a symbolic token (fs, del, dup, ins, *, =),
## or NA when the string carries no alt (e.g. "A298").
.acmg_parse_aa_change <- function(aa_change) {
  empty <- list(ref = NA_character_, pos = NA_integer_, alt = NA_character_)
  if (is.null(aa_change) || length(aa_change) != 1L || is.na(aa_change)) {
    return(empty)
  }
  s <- sub("^p\\.", "", as.character(aa_change))
  m <- regmatches(s, regexec("^([A-Za-z\\*])([0-9]+)([A-Za-z\\*=]*)$", s))[[1]]
  if (length(m) == 0L) {
    ## Fall back to a bare position (e.g. "175").
    mp <- regmatches(s, regexec("([0-9]+)", s))[[1]]
    if (length(mp) >= 1L) {
      return(list(ref = NA_character_,
                  pos = as.integer(mp[1]),
                  alt = NA_character_))
    }
    return(empty)
  }
  alt <- m[4]
  list(ref = toupper(m[2]),
       pos = as.integer(m[3]),
       alt = if (nzchar(alt)) alt else NA_character_)
}
