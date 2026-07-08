## Tests for the local data manifest: available_genes() and sha256
## checksum verification. These exercise only local files -- no network.

# Build a temporary cache dir populated with the two fixture bundles
# and a MANIFEST.tsv. `sha` controls what checksums the manifest lists:
#   "correct" -> real sha256 of each bundle
#   "wrong"   -> a deliberately bogus sha256 for TP53 (WDR31 stays real)
# Returns the versioned cache subdir (where bundles + manifest live).
.setup_cache_with_manifest <- function(sha = c("correct", "wrong")) {
  sha <- match.arg(sha)
  tmp <- tempfile()
  d <- file.path(tmp, "0.1.0")
  dir.create(d, recursive = TRUE)
  file.copy(test_path("fixtures", "TP53.rds"),  file.path(d, "TP53.rds"))
  file.copy(test_path("fixtures", "WDR31.rds"), file.path(d, "WDR31.rds"))

  tp53_sha  <- .sha256_file(file.path(d, "TP53.rds"))
  wdr31_sha <- .sha256_file(file.path(d, "WDR31.rds"))
  if (sha == "wrong") {
    tp53_sha <- paste(rep("0", nchar(tp53_sha)), collapse = "")
  }

  manifest <- data.frame(
    file   = c("TP53.rds", "WDR31.rds"),
    size   = file.info(file.path(d, c("TP53.rds", "WDR31.rds")))$size,
    sha256 = c(tp53_sha, wdr31_sha),
    stringsAsFactors = FALSE
  )
  write.table(manifest, file.path(d, "MANIFEST.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  Sys.setenv(MSAVARIANT_CACHE = tmp)
  d
}

test_that("available_genes() lists the genes in the manifest (TP53, WDR31)", {
  old <- Sys.getenv("MSAVARIANT_CACHE")
  on.exit(Sys.setenv(MSAVARIANT_CACHE = old))
  .setup_cache_with_manifest("correct")

  g <- available_genes()
  expect_type(g, "character")
  expect_setequal(g, c("TP53", "WDR31"))
})

test_that("available_genes() returns empty vector when no manifest present", {
  old <- Sys.getenv("MSAVARIANT_CACHE")
  on.exit(Sys.setenv(MSAVARIANT_CACHE = old))
  Sys.setenv(MSAVARIANT_CACHE = tempfile())   # empty, no manifest
  expect_identical(available_genes(), character(0))
})

test_that("checksum verification passes for a correct bundle", {
  old <- Sys.getenv("MSAVARIANT_CACHE")
  on.exit(Sys.setenv(MSAVARIANT_CACHE = old))
  d <- .setup_cache_with_manifest("correct")

  expect_true(.verify_checksum("TP53",  file.path(d, "TP53.rds")))
  expect_true(.verify_checksum("WDR31", file.path(d, "WDR31.rds")))
})

test_that("checksum verification fails for a corrupted bundle", {
  old <- Sys.getenv("MSAVARIANT_CACHE")
  on.exit(Sys.setenv(MSAVARIANT_CACHE = old))
  d <- .setup_cache_with_manifest("correct")

  # Corrupt the cached TP53 bundle by appending bytes; its sha256 now
  # differs from the (correct) manifest entry.
  con <- file(file.path(d, "TP53.rds"), open = "ab")
  writeBin(as.raw(c(0xDE, 0xAD, 0xBE, 0xEF)), con)
  close(con)

  expect_false(.verify_checksum("TP53", file.path(d, "TP53.rds")))
})

test_that("checksum verification is skipped (TRUE) when no manifest present", {
  old <- Sys.getenv("MSAVARIANT_CACHE")
  on.exit(Sys.setenv(MSAVARIANT_CACHE = old))
  tmp <- tempfile()
  d <- file.path(tmp, "0.1.0"); dir.create(d, recursive = TRUE)
  file.copy(test_path("fixtures", "TP53.rds"), file.path(d, "TP53.rds"))
  Sys.setenv(MSAVARIANT_CACHE = tmp)   # no MANIFEST.tsv written

  # No manifest -> nothing to check against -> permissive TRUE.
  expect_true(.verify_checksum("TP53", file.path(d, "TP53.rds")))
})

test_that("fetch_gene_data invalidates a cached bundle on checksum mismatch", {
  old <- Sys.getenv("MSAVARIANT_CACHE")
  on.exit(Sys.setenv(MSAVARIANT_CACHE = old))
  d <- .setup_cache_with_manifest("wrong")   # manifest lists bad TP53 sha
  cache_file <- file.path(d, "TP53.rds")
  expect_true(file.exists(cache_file))

  # The cached file is schema-valid but its sha256 disagrees with the
  # manifest -> fetch should warn, remove it, and return NULL.
  expect_warning(out <- fetch_gene_data("TP53"), "checksum")
  expect_null(out)
  expect_false(file.exists(cache_file))

  # WDR31's manifest sha is correct, so it loads fine from the same cache.
  expect_type(fetch_gene_data("WDR31"), "list")
})
