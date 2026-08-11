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

MSAVARIANT_DATA_VERSION <- "0.1.0"
MSAVARIANT_ZENODO_RECORD <- "PENDING_RECORD_ID"
MSAVARIANT_DATA_DOI <- "PENDING_DOI"

## ----- URL pattern -----------------------------------------------
## Host defaults to production zenodo.org. Set the MSAVARIANT_ZENODO_HOST
## env var (e.g. "https://sandbox.zenodo.org") to point fetches at the
## disposable sandbox for testing. Trailing slashes are tolerated.
.zenodo_host <- function() {
    sub("/+$", "", Sys.getenv("MSAVARIANT_ZENODO_HOST", "https://zenodo.org"))
}

## `name` is the file basename to fetch (a gene in per-gene mode, or a
## group in grouped mode); the deposit stores <name>.rds either way.
.zenodo_url <- function(name) {
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
        .zenodo_host(), MSAVARIANT_ZENODO_RECORD, name
    )
}

## Fetch a <name>.rds data file into `dest` (a group file in grouped
## mode, a gene file in per-gene mode). Returns TRUE on success, FALSE
## on failure. Source is the Zenodo deposit, unless MSAVARIANT_LOCAL_SOURCE
## points at a directory holding <name>.rds files -- an offline override
## used for air-gapped use and for testing the fetch path without network.
.download_data_file <- function(name, dest, quiet = FALSE) {
    tmp <- tempfile(fileext = ".rds")
    on.exit(unlink(tmp), add = TRUE)

    local_src <- Sys.getenv("MSAVARIANT_LOCAL_SOURCE", "")
    if (nzchar(local_src)) {
        src <- file.path(local_src, paste0(name, ".rds"))
        if (!file.exists(src)) {
            rlang::warn(sprintf(
                "MSAVARIANT_LOCAL_SOURCE set but '%s' not found.", src
            ))
            return(FALSE)
        }
        if (!file.copy(src, tmp, overwrite = TRUE)) {
            return(FALSE)
        }
    } else {
        url <- .zenodo_url(name)
        result <- tryCatch(
            utils::download.file(url, destfile = tmp, mode = "wb", quiet = TRUE),
            error = function(e) {
                rlang::warn(c(
                    sprintf("Could not download data file '%s'.", name),
                    "x" = conditionMessage(e),
                    "i" = "Check your internet connection."
                ))
                NULL
            }
        )
        if (is.null(result) || !file.exists(tmp) || file.info(tmp)$size == 0L) {
            return(FALSE)
        }
    }

    if (!quiet) {
        size_mb <- file.info(tmp)$size / 1024 / 1024
        if (size_mb > 2) message(sprintf("  downloaded %.1f MB.", size_mb))
    }
    ok <- file.rename(tmp, dest)
    if (!ok) file.copy(tmp, dest, overwrite = TRUE)
    TRUE
}

## ----- Cache location --------------------------------------------
## Platform-appropriate:
##   Linux:   ~/.cache/msaVariant/
##   macOS:   ~/Library/Caches/msaVariant/
##   Windows: %LOCALAPPDATA%/msaVariant/Cache/
## Override with MSAVARIANT_CACHE env var.

.cache_dir <- function() {
    env <- Sys.getenv("MSAVARIANT_CACHE", "")
    if (nzchar(env)) {
        return(env)
    }
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
#' Downloads (or reads from local cache) the annotation bundle for a
#' gene from the Zenodo deposit. Returns the deserialized list as
#' described in `DATA_FORMAT_SPEC.md`.
#'
#' Most users will not call this directly; the `get_domains()`,
#' `get_clinvar()`, etc. helpers and the `geom_*()` layers route
#' through it transparently.
#'
#' The public interface is the same regardless of how the deposit is
#' organised. When a gene->group index ships in the package
#' (`inst/extdata/gene_group_index.tsv`), the deposit packs many gene
#' bundles into a smaller number of "group" files; this function then
#' downloads the whole group once, caches it, and extracts the one
#' gene's bundle — so a second gene in the same group needs no further
#' download. Set the env var `MSAVARIANT_GROUPED` to a false-y value to
#' force the legacy one-file-per-gene mode. A bundle already present in
#' the cache under `<gene>.rds` (e.g. from [import_local_bundle()])
#' always takes precedence over grouped fetching.
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
#' @examples
#' ## Runnable with the shipped synthetic DEMO1 bundle (no network).
#' ## A temporary cache keeps the example off your real cache directory.
#' Sys.setenv(MSAVARIANT_CACHE = tempfile("msaVariant_cache_"))
#' import_local_bundle(
#'     system.file("extdata", "DEMO1.rds", package = "msaVariant"),
#'     gene = "DEMO1"
#' )
#' b <- fetch_gene_data("DEMO1")
#' names(b)
#'
#' \dontrun{
#' ## Real genes are downloaded from the Zenodo data deposit.
#' b <- fetch_gene_data("TP53")
#' }
#' @export
fetch_gene_data <- function(gene,
                            force_refresh = FALSE,
                            validate = TRUE,
                            verify_checksum = TRUE,
                            quiet = FALSE) {
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

    # A per-gene cache file takes precedence over grouped mode. This is
    # what import_local_bundle() writes, so an explicitly-imported local
    # bundle always wins and the legacy per-gene flow is preserved.
    if (!force_refresh && file.exists(cache_file)) {
        out <- .safe_read_rds(cache_file)
        if (!is.null(out) && validate) .validate_or_invalidate(out, cache_file)
        if (!is.null(out) && verify_checksum) {
            out <- .checksum_or_invalidate(gene, cache_file, out)
        }
        return(out)
    }

    # Grouped mode: if enabled and the gene is in the group index, fetch
    # the whole group once and extract this gene's bundle from it.
    group <- if (.grouped_enabled()) .gene_group(gene) else NA_character_
    if (!is.na(group)) {
        return(.fetch_grouped(
            gene, group, force_refresh, validate, verify_checksum, quiet
        ))
    }

    # Fallback: legacy per-gene download (gene not in any group).
    if (!quiet) {
        message(sprintf("msaVariant: downloading annotation for %s ...", gene))
    }
    if (!.download_data_file(gene, cache_file, quiet)) {
        rlang::warn(c(
            sprintf("Could not fetch annotation for %s.", gene),
            "i" = "If you have the data locally, pass it via the `data` argument of each geom."
        ))
        return(NULL)
    }

    out <- .safe_read_rds(cache_file)
    if (!is.null(out) && validate) .validate_or_invalidate(out, cache_file)
    if (!is.null(out) && verify_checksum) {
        out <- .checksum_or_invalidate(gene, cache_file, out)
    }
    out
}

# Grouped fetch: resolve the group file (download once, cache), then
# extract the requested gene's bundle from the group's named list. A
# second gene in the same group is served entirely from the cached group
# file with no further download. Returns the 7-element bundle, or NULL.
.fetch_grouped <- function(gene, group,
                           force_refresh = FALSE,
                           validate = TRUE,
                           verify_checksum = TRUE,
                           quiet = FALSE) {
    group_file <- .group_cache_path(group)

    # Drop a cached group that fails its manifest checksum, so it will be
    # re-fetched below.
    if (!force_refresh && file.exists(group_file) && verify_checksum) {
        if (!.verify_checksum(group, group_file)) {
            rlang::warn(c(
                sprintf("Cached group %s failed checksum verification; removing.", group),
                "x" = "sha256 does not match the local MANIFEST.tsv entry."
            ))
            unlink(group_file)
        }
    }

    if (force_refresh || !file.exists(group_file)) {
        if (!quiet) {
            message(sprintf(
                "msaVariant: downloading group %s (for %s) ...", group, gene
            ))
        }
        if (!.download_data_file(group, group_file, quiet)) {
            rlang::warn(sprintf("Could not fetch group '%s' for gene '%s'.", group, gene))
            return(NULL)
        }
        if (verify_checksum && !.verify_checksum(group, group_file)) {
            rlang::warn(c(
                sprintf("Downloaded group %s failed checksum verification; removing.", group),
                "x" = "sha256 does not match the local MANIFEST.tsv entry."
            ))
            unlink(group_file)
            return(NULL)
        }
    }

    group_list <- .safe_read_rds(group_file)
    if (is.null(group_list)) {
        return(NULL)
    }
    if (!is.list(group_list) || is.null(names(group_list)) ||
        !gene %in% names(group_list)) {
        rlang::warn(sprintf(
            "Gene '%s' not found inside group '%s' (index/deposit mismatch).",
            gene, group
        ))
        return(NULL)
    }

    out <- group_list[[gene]]
    if (!is.null(out) && validate) {
        res <- validate_gene_data(out)
        if (!isTRUE(res$valid)) {
            rlang::warn(c(
                sprintf("Bundle for %s (from group %s) failed validation.", gene, group),
                "x" = paste(res$issues, collapse = "; ")
            ))
            return(NULL)
        }
    }
    out
}

# Verify a cached file against the local MANIFEST.tsv checksum. If a
# manifest checksum exists and does not match, remove the file and
# return NULL (next call will redownload). Returns `out` unchanged
# when there is no manifest / no entry / no checksum backend, keeping
# behaviour identical when no manifest is present.
.checksum_or_invalidate <- function(gene, cache_file, out) {
    if (!file.exists(cache_file)) {
        return(out)
    } # already invalidated
    if (.verify_checksum(gene, cache_file)) {
        return(out)
    }
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
        }
    )
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
#' @examples
#' Sys.setenv(MSAVARIANT_CACHE = tempfile("msaVariant_cache_"))
#' import_local_bundle(
#'     system.file("extdata", "DEMO1.rds", package = "msaVariant"),
#'     gene = "DEMO1"
#' )
#' clear_cache("DEMO1") # remove one gene
#' clear_cache() # remove everything
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
#' @examples
#' Sys.setenv(MSAVARIANT_CACHE = tempfile("msaVariant_cache_"))
#' cache_location()
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
#' @examples
#' Sys.setenv(MSAVARIANT_CACHE = tempfile("msaVariant_cache_"))
#' import_local_bundle(
#'     system.file("extdata", "DEMO1.rds", package = "msaVariant"),
#'     gene = "DEMO1"
#' )
#' cache_summary()
#' @export
cache_summary <- function() {
    base <- file.path(.cache_dir(), MSAVARIANT_DATA_VERSION)
    if (!dir.exists(base)) {
        return(data.frame(
            gene = character(0),
            size_kb = numeric(0),
            cached_on = as.Date(character(0)),
            stringsAsFactors = FALSE
        ))
    }
    files <- list.files(base, pattern = "\\.rds$", full.names = TRUE)
    if (length(files) == 0L) {
        return(data.frame(
            gene = character(0),
            size_kb = numeric(0),
            cached_on = as.Date(character(0)),
            stringsAsFactors = FALSE
        ))
    }
    info <- file.info(files)
    out <- data.frame(
        gene = sub("\\.rds$", "", basename(files)),
        size_kb = round(info$size / 1024, 1),
        cached_on = as.Date(info$mtime),
        stringsAsFactors = FALSE
    )
    out[order(-out$size_kb), , drop = FALSE]
}
