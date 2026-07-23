## -----------------------------------------------------------------
## msaVariant: Zenodo-backed per-gene fetcher (v0.1)
## -----------------------------------------------------------------
##
## ARCHITECTURE
## ------------
## One combined `.rds` per gene, fetched from Zenodo, cached locally.
## The file holds a named list of seven elements (meta, domains,
## clinvar, gnomad, alphamissense, revel, cadd) -- see DATA_FORMAT_SPEC.md
## for the authoritative schema.
##
## User-facing `get_*()` functions (in fetchers.R) extract the
## relevant slice from this combined object. The geoms in
## geom_annotations.R route through `get_*()` when given a gene
## name, or accept user-supplied data frames directly.

## ----- Configuration ---------------------------------------------
## When the data deposit is published, update these constants.
## Zenodo records are immutable per version, so pinning to a specific
## record ID is the reproducibility contract.

MSAVARIANT_DATA_VERSION  <- "0.1.0"
MSAVARIANT_ZENODO_RECORD <- "PENDING_RECORD_ID"
MSAVARIANT_DATA_DOI      <- "PENDING_DOI"

## ----- URL pattern -----------------------------------------------
## Host defaults to production zenodo.org. Set the MSAVARIANT_ZENODO_HOST
## env var (e.g. "https://sandbox.zenodo.org") to point fetches at the
## disposable sandbox for testing. Trailing slashes are tolerated.
.zenodo_host <- function() {
  sub("/+$", "", Sys.getenv("MSAVARIANT_ZENODO_HOST", "https://zenodo.org"))
}

.zenodo_url <- function(gene) {
  if (MSAVARIANT_ZENODO_RECORD == "PENDING_RECORD_ID") {
    rlang::abort(c(
      "msaVariant data deposit URL has not been configured.",
      "i" = "Edit R/fetch_core.R and set MSAVARIANT_ZENODO_RECORD",
      "i" = "and MSAVARIANT_DATA_DOI after uploading the data to Zenodo.",
      "i" = "See data-raw/ZENODO_UPLOAD.md for instructions."
    ))
  }
  ## Modern Zenodo (InvenioRDM) serves files under the plural /records/
  ## path; the legacy /record/ singular still 301-redirects to it.
  sprintf(
    "%s/records/%s/files/%s.rds?download=1",
    .zenodo_host(), MSAVARIANT_ZENODO_RECORD, gene
  )
}

## ----- Cache location --------------------------------------------
## Platform-appropriate:
##   Linux:   ~/.cache/msaVariant/
##   macOS:   ~/Library/Caches/msaVariant/
##   Windows: %LOCALAPPDATA%/msaVariant/Cache/
## Override with MSAVARIANT_CACHE env var.

.cache_dir <- function() {
  env <- Sys.getenv("MSAVARIANT_CACHE", "")
  if (nzchar(env)) return(env)
  tools::R_user_dir("msaVariant", which = "cache")
}

.cache_path <- function(gene) {
  d <- file.path(.cache_dir(), MSAVARIANT_DATA_VERSION)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  file.path(d, paste0(gene, ".rds"))
}

## ----- The main fetch function -----------------------------------

#' Fetch the combined annotation file for a gene
#'
#' Downloads (or reads from local cache) the per-gene annotation
#' bundle from the Zenodo deposit. Returns the deserialized list
#' as described in `DATA_FORMAT_SPEC.md`.
#'
#' Most users will not call this directly; the `get_domains()`,
#' `get_clinvar()`, etc. helpers and the `geom_*()` layers route
#' through it transparently.
#'
#' @param gene HGNC gene symbol (e.g. `"PATL1"`).
#' @param force_refresh If `TRUE`, redownload even if cached.
#' @param validate If `TRUE` (default), validate the file against
#'   the package's format spec before returning.
#' @param verify_checksum If `TRUE` (default), and a local
#'   `MANIFEST.tsv` is present, verify the bundle's sha256 against the
#'   manifest entry; a mismatched file is treated as corrupt, removed
#'   from the cache, and `NULL` is returned. Has no effect when no
#'   manifest is present.
#' @param quiet If `TRUE`, suppress "Downloading..." messages.
#' @return A named list with 7 elements (`meta`, `domains`,
#'   `clinvar`, `gnomad`, `alphamissense`, `revel`, `cadd`), or
#'   `NULL` with a warning on failure.
#' @export
fetch_gene_data <- function(gene,
                             force_refresh   = FALSE,
                             validate        = TRUE,
                             verify_checksum = TRUE,
                             quiet           = FALSE) {
  if (!is.character(gene) || length(gene) != 1L || !nzchar(gene)) {
    rlang::abort("`gene` must be a single non-empty HGNC symbol.")
  }
  # Reject suspicious symbols. HGNC symbols are alphanumeric plus
  # a few permitted punctuation marks; reject anything that looks
  # like a path-traversal attempt.
  if (!grepl("^[A-Za-z0-9][A-Za-z0-9._-]*$", gene)) {
    rlang::abort(sprintf("Suspicious gene symbol '%s' rejected.", gene))
  }

  cache_file <- .cache_path(gene)

  if (!force_refresh && file.exists(cache_file)) {
    out <- .safe_read_rds(cache_file)
    if (!is.null(out) && validate) .validate_or_invalidate(out, cache_file)
    if (!is.null(out) && verify_checksum)
      out <- .checksum_or_invalidate(gene, cache_file, out)
    return(out)
  }

  url <- .zenodo_url(gene)

  # Friendly progress message scaled to expected size. We don't
  # know the actual size until the HEAD request lands, but for a
  # generic message we just say "downloading".
  if (!quiet) {
    message(sprintf("msaVariant: downloading annotation for %s ...", gene))
  }

  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)

  result <- tryCatch(
    utils::download.file(url, destfile = tmp, mode = "wb", quiet = TRUE),
    error = function(e) {
      rlang::warn(c(
        sprintf("Could not download annotation for %s.", gene),
        "x" = conditionMessage(e),
        "i" = "Check your internet connection.",
        "i" = "If you have the data locally, pass it via the `data` argument of each geom."
      ))
      NULL
    }
  )
  if (is.null(result) || !file.exists(tmp) || file.info(tmp)$size == 0L) {
    return(NULL)
  }

  # Friendly size message after download
  if (!quiet) {
    size_mb <- file.info(tmp)$size / 1024 / 1024
    if (size_mb > 2) {
      message(sprintf("  downloaded %.1f MB.", size_mb))
    }
  }

  # Move from temp to cache (atomic on same filesystem)
  ok <- file.rename(tmp, cache_file)
  if (!ok) file.copy(tmp, cache_file, overwrite = TRUE)

  out <- .safe_read_rds(cache_file)
  if (!is.null(out) && validate) .validate_or_invalidate(out, cache_file)
  if (!is.null(out) && verify_checksum)
    out <- .checksum_or_invalidate(gene, cache_file, out)
  out
}

# Verify a cached file against the local MANIFEST.tsv checksum. If a
# manifest checksum exists and does not match, remove the file and
# return NULL (next call will redownload). Returns `out` unchanged
# when there is no manifest / no entry / no checksum backend, keeping
# behaviour identical when no manifest is present.
.checksum_or_invalidate <- function(gene, cache_file, out) {
  if (!file.exists(cache_file)) return(out)   # already invalidated
  if (.verify_checksum(gene, cache_file)) return(out)
  rlang::warn(c(
    sprintf("Cached %s bundle failed checksum verification; removing.", gene),
    "x" = "sha256 does not match the local MANIFEST.tsv entry.",
    "i" = "The file may be corrupt; it will be re-fetched on next call."
  ))
  unlink(cache_file)
  NULL
}

# Validate a cached file. If it fails, remove it and return NULL --
# next call will redownload.
.validate_or_invalidate <- function(obj, cache_file) {
  res <- validate_gene_data(obj)
  if (!isTRUE(res$valid)) {
    rlang::warn(c(
      sprintf("Cached gene file at %s failed validation; removing.", cache_file),
      "x" = paste(res$issues, collapse = "; ")
    ))
    unlink(cache_file)
    return(invisible(FALSE))
  }
  invisible(TRUE)
}

# Read RDS with a clear error if the file is corrupt.
.safe_read_rds <- function(path) {
  tryCatch(readRDS(path),
           error = function(e) {
             rlang::warn(c(
               sprintf("Cached file at %s is unreadable; removing.", path),
               "x" = conditionMessage(e)
             ))
             unlink(path)
             NULL
           })
}

# Null-coalescing helper used above.
`%||%` <- function(x, y) if (is.null(x)) y else x

## ----- Cache management ------------------------------------------

#' Clear the msaVariant download cache
#'
#' Removes locally-cached annotation files. Use this if you suspect
#' a cached file is stale or corrupt, or if you want to free disk
#' space.
#'
#' @param gene If supplied, only that gene's cached file is removed.
#'   If `NULL` (default), the entire cache is removed.
#' @return Invisibly, the number of files deleted.
#' @export
clear_cache <- function(gene = NULL) {
  base <- file.path(.cache_dir(), MSAVARIANT_DATA_VERSION)
  if (!dir.exists(base)) {
    message("Cache is already empty.")
    return(invisible(0L))
  }
  if (!is.null(gene)) {
    f <- .cache_path(gene)
    if (!file.exists(f)) {
      message(sprintf("No cache file for '%s'.", gene))
      return(invisible(0L))
    }
    unlink(f)
    message(sprintf("Cleared cached file for %s.", gene))
    return(invisible(1L))
  }
  files <- list.files(base, full.names = TRUE, recursive = TRUE)
  unlink(base, recursive = TRUE)
  message(sprintf("Cleared %d cached file(s).", length(files)))
  invisible(length(files))
}

#' Show the location of the msaVariant cache
#' @return The cache directory path (character).
#' @export
cache_location <- function() {
  d <- .cache_dir()
  if (!dir.exists(d)) {
    message("Cache directory does not exist yet (no annotations have been downloaded).")
  }
  d
}

#' Summarize what's currently cached
#'
#' Lists all locally-cached per-gene annotation files, with their
#' size and date of caching. Useful for inspecting disk use after
#' querying many genes.
#'
#' @return A `data.frame` with columns `gene`, `size_kb`, `cached_on`,
#'   sorted by size descending. Returns an empty data.frame if the
#'   cache is empty.
#' @export
cache_summary <- function() {
  base <- file.path(.cache_dir(), MSAVARIANT_DATA_VERSION)
  if (!dir.exists(base)) {
    return(data.frame(gene = character(0),
                      size_kb = numeric(0),
                      cached_on = as.Date(character(0)),
                      stringsAsFactors = FALSE))
  }
  files <- list.files(base, pattern = "\\.rds$", full.names = TRUE)
  if (length(files) == 0L) {
    return(data.frame(gene = character(0),
                      size_kb = numeric(0),
                      cached_on = as.Date(character(0)),
                      stringsAsFactors = FALSE))
  }
  info <- file.info(files)
  out <- data.frame(
    gene      = sub("\\.rds$", "", basename(files)),
    size_kb   = round(info$size / 1024, 1),
    cached_on = as.Date(info$mtime),
    stringsAsFactors = FALSE
  )
  out[order(-out$size_kb), , drop = FALSE]
}
