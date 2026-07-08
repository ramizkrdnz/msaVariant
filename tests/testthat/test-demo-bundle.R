## Smoke tests for the shipped synthetic DEMO1 example bundle.
## These guard the out-of-the-box Quick Example / man example.

test_that("shipped DEMO1 bundle is present, valid, and clearly synthetic", {
  bundle_path <- system.file("extdata", "DEMO1.rds", package = "msaVariant")
  expect_true(nzchar(bundle_path))
  expect_true(file.exists(bundle_path))

  b <- readRDS(bundle_path)
  v <- validate_gene_data(b, strict = FALSE)
  expect_true(isTRUE(v$valid))

  # Labelled synthetic so nobody mistakes it for real predictions
  expect_identical(as.character(b$meta$gene), "DEMO1")
  expect_match(b$meta$source_versions, "SYNTHETIC", ignore.case = TRUE)

  # Every annotation track is populated (so the demo figure looks full)
  for (tbl in c("domains", "clinvar", "gnomad",
                "alphamissense", "revel", "cadd")) {
    expect_gt(nrow(b[[tbl]]), 0L)
  }
})

test_that("DEMO1 p.R21H fires all five ACMG codes", {
  b <- readRDS(system.file("extdata", "DEMO1.rds", package = "msaVariant"))
  codes <- compute_acmg_codes(b, 21, "R21H")
  expect_identical(codes, c("PS1", "PM1", "PM2", "PM5", "PP3"))
})

test_that("DEMO1 renders end-to-end via plot_variant_overlay", {
  withr::local_envvar(MSAVARIANT_CACHE = tempfile("msaVariant_cache_"))
  import_local_bundle(
    system.file("extdata", "DEMO1.rds", package = "msaVariant"),
    gene = "DEMO1"
  )
  fa <- system.file("extdata", "demo_aligned.fasta", package = "msaVariant")
  p <- plot_variant_overlay(
    gene = "DEMO1", aligned_fasta = fa,
    variant_pos = 21, variant_label = "p.R21H"
  )
  expect_s3_class(p, "patchwork")
})
