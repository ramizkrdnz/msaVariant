## Layout-planning tests for plot_variant_overlay().
## We assert on the assembled layout structure (number of active
## label/data rows in the returned patchwork) and on the strict
## separation between ACMG computation and ACMG *display* — not on
## pixel output. Building the plot object (no ggsave) is enough to
## exercise the row-planning logic.

## The plot reads the bundle from the package cache, so seed it from the
## fixture; the aligned FASTA fixture is passed explicitly.
import_local_bundle(test_path("fixtures", "TP53.rds"))
fasta <- test_path("fixtures", "tp53_aligned.fasta")

## Build a figure with the standard TP53 R175H inputs plus any overrides.
build <- function(...) {
  suppressWarnings(
    plot_variant_overlay(gene = "TP53", aligned_fasta = fasta,
                         variant_pos = 175, variant_label = "p.R175H",
                         window = c(160, 200), ...)
  )
}

## A patchwork built by wrap_plots stores the first plot as the base and
## the remaining plots in $patches$plots; each layout row contributes a
## (label, data) pair, so rows = panels / 2.
n_panels <- function(p) length(p$patches$plots) + 1L
n_rows   <- function(p) n_panels(p) / 2

## Full stack: variant + ACMG + MSA + 6 tracks + legend dashboard = 10 rows.
FULL_ROWS <- 10L

test_that("all layers ON assembles the full expected number of rows", {
  p <- build()
  expect_s3_class(p, "patchwork")
  expect_equal(n_panels(p) %% 2L, 0L)       # always paired (label,data)
  expect_equal(n_rows(p), FULL_ROWS)
})

test_that("turning off each layer removes exactly one row (no error)", {
  full <- n_rows(build())
  toggles <- c("show_domain", "show_clinvar", "show_gnomad",
               "show_cadd", "show_alphamissense", "show_revel", "show_acmg")
  for (tg in toggles) {
    args <- stats::setNames(list(FALSE), tg)
    p <- do.call(build, args)
    expect_equal(n_rows(p), full - 1L,
                 info = sprintf("%s=FALSE should drop exactly one row", tg))
    expect_s3_class(p, "patchwork")
  }
})

test_that("minimal case (MSA + variant + ACMG only) assembles cleanly", {
  p <- build(show_domain = FALSE, show_clinvar = FALSE, show_gnomad = FALSE,
             show_cadd = FALSE, show_alphamissense = FALSE, show_revel = FALSE)
  ## Exactly three rows, no orphan/blank rows (panels are all paired).
  expect_equal(n_rows(p), 3L)
  expect_equal(n_panels(p) %% 2L, 0L)
  expect_s3_class(p, "patchwork")
})

test_that("turning off two layers at once removes exactly two rows", {
  full <- n_rows(build())
  p <- build(show_gnomad = FALSE, show_cadd = FALSE)
  expect_equal(n_rows(p), full - 2L)
})

test_that("acmg_codes filtering never changes the computed codes", {
  b <- fetch_gene_data("TP53")
  all_codes <- compute_acmg_codes(b, 175, "R175H")
  expect_setequal(all_codes, c("PS1", "PM1", "PM2", "PM5", "PP3"))

  ## Whatever we ask to *display*, the underlying computation is identical.
  for (subset in list(c("PS1"), c("PS1", "PM1"), c("PP3"), character(0),
                      c("PS1","PM1","PM2","PM5","PP3"))) {
    p <- build(acmg_codes = subset)          # display filter only
    expect_s3_class(p, "patchwork")
    ## The ACMG row is still present regardless of how many chips show.
    expect_equal(n_rows(p), 10L)
    ## Recompute -> unchanged, proving display != logic.
    expect_setequal(compute_acmg_codes(b, 175, "R175H"), all_codes)
  }
})
