## -----------------------------------------------------------------
## msaVariant: local data manifest (MANIFEST.tsv) support
## -----------------------------------------------------------------
##
## The data deposit ships a `MANIFEST.tsv` alongside the per-gene
## `.rds` bundles (columns: file, size, sha256; see
## data-raw/build_scripts/build_manifest.R). At runtime the manifest
## lives next to the cached bundles, at
##   <cache>/<data version>/MANIFEST.tsv
## and gives us two things the bare URL scheme could not:
##   (1) the set of genes actually available, and
##   (2) an integrity check (sha256) for a fetched/imported bundle.
##
## Everything here degrades gracefully: with no manifest present the
## helpers are no-ops and package behaviour is exactly as before.

## ----- Manifest location & reading -------------------------------

# Path to the local manifest, alongside the cached bundles.
.manifest_path <- function() {
  file.path(.cache_dir(), MSAVARIANT_DATA_VERSION, "MANIFEST.tsv")
}

# Read the manifest as a data.frame, or NULL if absent/unreadable.
# Expected columns: file, size, sha256.
.read_manifest <- function(path = .manifest_path()) {
  if (!file.exists(path)) return(NULL)
  m <- tryCatch(
    utils::read.delim(path, sep = "\t", stringsAsFactors = FALSE,
                      colClasses = "character"),
    error = function(e) {
      rlang::warn(c(
        sprintf("MANIFEST at %s is unreadable; ignoring.", path),
        "x" = conditionMessage(e)
      ))
      NULL
    }
  )
  if (is.null(m) || !all(c("file", "sha256") %in% names(m))) {
    return(NULL)
  }
  m
}

# The expected sha256 for a gene from the manifest, or NA_character_
# if the manifest is absent or has no entry / no checksum for it.
.manifest_sha256 <- function(gene, manifest = .read_manifest()) {
  if (is.null(manifest)) return(NA_character_)
  row <- manifest[manifest$file == paste0(gene, ".rds"), , drop = FALSE]
  if (nrow(row) < 1L) return(NA_character_)
  sha <- row$sha256[1]
  if (is.null(sha) || is.na(sha) || !nzchar(sha)) return(NA_character_)
  sha
}

## ----- Checksum backend ------------------------------------------

# sha256 of a file. Uses the 'digest' package if available, else a
# system tool (sha256sum / shasum -a 256). Returns NA_character_ if no
# backend is available, so callers can treat verification as skipped.
.sha256_file <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  if (requireNamespace("digest", quietly = TRUE)) {
    return(tryCatch(digest::digest(file = path, algo = "sha256"),
                    error = function(e) NA_character_))
  }
  if (nzchar(Sys.which("sha256sum"))) {
    out <- tryCatch(system2("sha256sum", path, stdout = TRUE, stderr = FALSE),
                    error = function(e) NA_character_)
  } else if (nzchar(Sys.which("shasum"))) {
    out <- tryCatch(system2("shasum", c("-a", "256", path),
                            stdout = TRUE, stderr = FALSE),
                    error = function(e) NA_character_)
  } else {
    return(NA_character_)
  }
  if (length(out) < 1L || is.na(out[1])) return(NA_character_)
  sub("\\s.*$", "", out[1])
}

# Verify a bundle file against its manifest checksum.
# Returns TRUE  -> matches, OR verification not applicable (no manifest,
#                  no entry for this gene, or no checksum backend).
#         FALSE -> a manifest checksum exists and does NOT match.
# The permissive default on "not applicable" is what keeps behaviour
# unchanged when no manifest is present.
.verify_checksum <- function(gene, path, manifest = .read_manifest()) {
  expected <- .manifest_sha256(gene, manifest)
  if (is.na(expected)) return(TRUE)          # nothing to check against
  actual <- .sha256_file(path)
  if (is.na(actual)) return(TRUE)            # no backend -> skip
  identical(tolower(actual), tolower(expected))
}

## ----- Public API ------------------------------------------------

#' List genes available in the local data manifest
#'
#' Reads the `MANIFEST.tsv` that ships alongside the cached per-gene
#' bundles (see [cache_location()]) and returns the gene symbols it
#' lists. This is the set of genes for which annotation data can be
#' loaded without a further download.
#'
#' If no manifest is present (for example before any data has been
#' downloaded or imported), an empty character vector is returned.
#'
#' @return A sorted character vector of HGNC gene symbols. Empty if
#'   no manifest is found.
#' @seealso [fetch_gene_data()], [import_local_bundle()],
#'   [cache_summary()]
#' @export
#' @examples
#' \dontrun{
#' available_genes()
#' }
available_genes <- function() {
  m <- .read_manifest()
  if (is.null(m)) return(character(0))
  genes <- sub("\\.rds$", "", m$file)
  genes <- genes[nzchar(genes)]
  sort(unique(genes))
}
