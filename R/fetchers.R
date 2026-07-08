## -----------------------------------------------------------------
## msaVariant: public per-annotation fetchers (v0.1)
## -----------------------------------------------------------------
##
## Each `get_*("GENE")` function returns the matching slice of the
## combined per-gene object fetched from Zenodo. Cached locally
## after first use.
##
## Returns NULL with a warning if the underlying download failed
## (so callers can fall back gracefully).

.fetch_slice <- function(gene, table_name, force_refresh = FALSE) {
  obj <- fetch_gene_data(gene, force_refresh = force_refresh)
  if (is.null(obj)) return(NULL)
  if (!table_name %in% names(obj)) {
    rlang::warn(sprintf("'%s' table missing from %s gene data.",
                        table_name, gene))
    return(NULL)
  }
  obj[[table_name]]
}

#' Get InterPro/Pfam domains for a gene
#'
#' Returns the `domains` slice of the per-gene Zenodo bundle.
#'
#' @param gene HGNC gene symbol (e.g. `"PATL1"`).
#' @param force_refresh Redownload even if cached.
#' @return A `data.frame` with columns `start`, `end`, `name`,
#'   `accession`, `source` (factor), and optionally `evidence`.
#'   Returns `NULL` if download failed.
#' @export
get_domains <- function(gene, force_refresh = FALSE) {
  .fetch_slice(gene, "domains", force_refresh)
}

#' Get ClinVar variants for a gene
#'
#' Returns the `clinvar` slice. See `DATA_FORMAT_SPEC.md` for column
#' definitions.
#'
#' @param gene HGNC gene symbol.
#' @param force_refresh Redownload even if cached.
#' @return A `data.frame` with `pos`, `aa_ref`, `aa_alt`,
#'   `aa_change`, `significance` (factor), `review_status` (factor),
#'   `clinvar_id`, and optionally `condition`, `last_evaluated`.
#'   Returns `NULL` if download failed.
#' @export
get_clinvar <- function(gene, force_refresh = FALSE) {
  .fetch_slice(gene, "clinvar", force_refresh)
}

#' Get gnomAD per-variant allele frequencies for a gene
#'
#' Returns the `gnomad` slice with one row per coding variant in
#' gnomAD v4.1 (or whatever version was current at deposit build
#' time; see `attr(., "source_versions")`).
#'
#' @param gene HGNC gene symbol.
#' @param force_refresh Redownload even if cached.
#' @return A `data.frame` per `DATA_FORMAT_SPEC.md`. Returns `NULL`
#'   if download failed.
#' @export
get_gnomad <- function(gene, force_refresh = FALSE) {
  .fetch_slice(gene, "gnomad", force_refresh)
}

#' Get AlphaMissense per-substitution scores for a gene
#'
#' Returns the full per-substitution table (every missense
#' substitution at every residue, ~19 rows per residue).
#'
#' Note: AlphaMissense is licensed CC-BY-NC-SA 4.0. Commercial
#' use is restricted.
#'
#' @param gene HGNC gene symbol.
#' @param force_refresh Redownload even if cached.
#' @return A `data.frame` with `pos`, `aa_ref`, `aa_alt`, `aa_change`,
#'   `am_score`, `am_class` (factor: likely_benign / ambiguous /
#'   likely_pathogenic). Returns `NULL` if download failed.
#' @export
get_alphamissense <- function(gene, force_refresh = FALSE) {
  .fetch_slice(gene, "alphamissense", force_refresh)
}

#' Get REVEL per-substitution scores for a gene
#'
#' Per-substitution REVEL ensemble scores in `[0, 1]`. Higher
#' indicates more likely pathogenic.
#'
#' @param gene HGNC gene symbol.
#' @param force_refresh Redownload even if cached.
#' @return A `data.frame` per spec; returns `NULL` on failure.
#' @export
get_revel <- function(gene, force_refresh = FALSE) {
  .fetch_slice(gene, "revel", force_refresh)
}

#' Get CADD per-substitution PHRED scores for a gene
#'
#' Per-substitution CADD PHRED-scaled deleteriousness scores.
#' Higher = more deleterious; PHRED 20 ~ top 1% deleterious.
#'
#' Note: CADD is licensed for **non-commercial use only**.
#'
#' @param gene HGNC gene symbol.
#' @param force_refresh Redownload even if cached.
#' @return A `data.frame` per spec; returns `NULL` on failure.
#' @export
get_cadd <- function(gene, force_refresh = FALSE) {
  .fetch_slice(gene, "cadd", force_refresh)
}

#' Get gene metadata
#'
#' Returns the `meta` slice — one row with UniProt ID, protein
#' length, Ensembl IDs, build date, and source-version stamps.
#'
#' @param gene HGNC gene symbol.
#' @param force_refresh Redownload even if cached.
#' @return A 1-row `data.frame`; returns `NULL` on failure.
#' @export
get_gene_meta <- function(gene, force_refresh = FALSE) {
  .fetch_slice(gene, "meta", force_refresh)
}
