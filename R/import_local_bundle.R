## ===================================================================
## msaVariant: import_local_bundle()
## ===================================================================
##
## Copy a locally-built `.rds` bundle into the package cache so that
## `fetch_gene_data()` and `plot_variant_overlay()` can find it.
## This is the convenient path for testing offline or before Zenodo
## hosting is wired up.

#' Import a locally-built gene bundle into the cache
#'
#' Copies the supplied `.rds` file into the package's cache directory
#' (see [cache_location()]) so that [fetch_gene_data()] and
#' [plot_variant_overlay()] can find it.
#'
#' @param path Path to a `.rds` file produced by your build script.
#'   The file must conform to the gene-data schema enforced by
#'   [validate_gene_data()].
#' @param gene Optional gene symbol; if `NULL` (default) the symbol
#'   is read from `bundle$meta$gene` inside the file. Used as the
#'   destination filename: `<GENE>.rds`.
#' @param overwrite Logical, default `TRUE`. If `FALSE` and a bundle
#'   for this gene is already present, an error is thrown.
#' @param verify_checksum Logical, default `TRUE`. If a local
#'   `MANIFEST.tsv` is present and lists this gene, the imported file's
#'   sha256 is compared against the manifest entry and a warning is
#'   emitted on mismatch. The file is still imported (you asked for it
#'   explicitly); the warning simply flags a possible corruption or a
#'   stale manifest. No effect when no manifest is present.
#'
#' @return Invisibly, the destination path of the imported bundle.
#'
#' @examples
#' \dontrun{
#' import_local_bundle("/path/to/TP53.rds")
#' import_local_bundle("/path/to/WDR31.rds")
#' p <- plot_variant_overlay(
#'   gene = "TP53",
#'   aligned_fasta = "tp53_aligned.fasta",
#'   variant_pos = 175
#' )
#' }
#'
#' @export
import_local_bundle <- function(path, gene = NULL, overwrite = TRUE,
                                 verify_checksum = TRUE) {
  if (!file.exists(path)) {
    stop("File not found: ", path)
  }
  bundle <- readRDS(path)
  
  ## Validate before importing — refuse non-spec files
  v <- validate_gene_data(bundle, strict = FALSE)
  if (!isTRUE(v$valid)) {
    stop("Bundle failed validation:\n  ",
         paste(v$issues, collapse = "\n  "))
  }
  
  if (is.null(gene)) {
    if (is.null(bundle$meta) || is.null(bundle$meta$gene)) {
      stop("Cannot infer gene symbol from bundle$meta$gene; pass `gene=`.")
    }
    gene <- as.character(bundle$meta$gene[1])
  }
  
  cache_root <- file.path(cache_location(), bundle$meta$data_version %||% "0.1.0")
  dir.create(cache_root, recursive = TRUE, showWarnings = FALSE)
  dest <- file.path(cache_root, paste0(gene, ".rds"))
  
  if (file.exists(dest) && !overwrite) {
    stop("Bundle already in cache: ", dest,
         " (set overwrite = TRUE to replace)")
  }
  file.copy(path, dest, overwrite = TRUE)

  ## Optional integrity check against a local MANIFEST.tsv (if any).
  if (isTRUE(verify_checksum) && !.verify_checksum(gene, dest)) {
    rlang::warn(c(
      sprintf("Imported %s bundle does not match its MANIFEST.tsv sha256.", gene),
      "x" = "The file's checksum differs from the manifest entry.",
      "i" = "The bundle was imported anyway; verify it is the intended build."
    ))
  }

  message("Imported ", gene, " bundle -> ", dest)
  invisible(dest)
}

## Internal infix fallback so the call above doesn't break on older R
`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1L && is.na(a))) b else a
