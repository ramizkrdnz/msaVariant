## Lock-in tests for compute_acmg_codes().
## These capture the behaviour already validated for the ACMG strip:
##   * TP53 R175H -> PS1, PM1, PM2, PM5, PP3 (PP3 fires at min = 1/2/3)
##   * WDR31 A298V -> PM1, PM2 only (PP3 predictors pass 0/3)
##   * CADD boundary: phred == 25 counts as passing (>= 25, not > 25)
## They must not require network access or the package cache; the two
## real bundles are read from tests/testthat/fixtures/.

tp53  <- readRDS(test_path("fixtures", "TP53.rds"))
wdr31 <- readRDS(test_path("fixtures", "WDR31.rds"))

test_that("TP53 R175H returns all five codes at pp3_min_predictors 1, 2, 3", {
  expected <- c("PS1", "PM1", "PM2", "PM5", "PP3")
  for (m in 1:3) {
    codes <- compute_acmg_codes(tp53, 175, "R175H", pp3_min_predictors = m)
    expect_setequal(codes, expected)
    expect_true("PP3" %in% codes,
                info = sprintf("PP3 should fire at pp3_min_predictors = %d", m))
    expect_length(codes, 5)
  }
})

test_that("TP53 R175H fires all five codes with the default (no pp3 arg)", {
  codes <- compute_acmg_codes(tp53, 175, "R175H")
  expect_setequal(codes, c("PS1", "PM1", "PM2", "PM5", "PP3"))
})

test_that("WDR31 A298V returns PM1 and PM2 but not PP3 at default", {
  codes <- compute_acmg_codes(wdr31, 298, "A298V")   # default pp3_min_predictors = 2
  expect_setequal(codes, c("PM1", "PM2"))
  expect_false("PP3" %in% codes)
})

## --- CADD boundary: phred exactly 25 must PASS (>= 25, not > 25) ------
## Build a minimal synthetic bundle where the ONLY PP3 predictor that
## can pass is CADD, so PP3 at pp3_min_predictors = 1 is a direct probe
## of the CADD comparison operator.
make_min_bundle <- function(cadd_phred) {
  empty <- function(cols) {
    df <- as.data.frame(setNames(rep(list(character(0)), length(cols)), cols))
    df
  }
  list(
    meta    = data.frame(gene = "TESTG", protein_length = 200L),
    domains = data.frame(start = integer(0), end = integer(0)),  # no PM1
    clinvar = data.frame(pos = integer(0), aa_alt = character(0),
                         aa_change = character(0), significance = character(0)),
    gnomad  = data.frame(aa_change = "A100V", af_joint = 0.5),   # common -> no PM2
    alphamissense = data.frame(aa_change = "A100V", am_class = "likely_benign"), # fails
    revel   = data.frame(aa_change = "A100V", revel_score = 0.50),               # fails
    cadd    = data.frame(aa_change = "A100V", cadd_phred = cadd_phred)
  )
}

test_that("CADD phred exactly 25 counts as passing (>= 25)", {
  b25 <- make_min_bundle(25)
  ## Only CADD can pass -> 1/3. At min = 1, PP3 fires iff CADD(25) passes.
  codes <- compute_acmg_codes(b25, 100, "A100V", pp3_min_predictors = 1)
  expect_true("PP3" %in% codes,
              info = "CADD phred == 25 must pass the >= 25 threshold")
})

test_that("CADD phred 24 does not pass (boundary is exactly 25)", {
  b24 <- make_min_bundle(24)
  codes <- compute_acmg_codes(b24, 100, "A100V", pp3_min_predictors = 1)
  expect_false("PP3" %in% codes,
               info = "CADD phred 24 is below 25 and must not pass")
})

test_that("With CADD == 25 as the sole passer, PP3 needs the >=25 count", {
  b25 <- make_min_bundle(25)
  ## 1 predictor passes; require 2 -> PP3 must NOT fire.
  expect_false("PP3" %in% compute_acmg_codes(b25, 100, "A100V", pp3_min_predictors = 2))
  ## require 1 -> PP3 fires.
  expect_true("PP3"  %in% compute_acmg_codes(b25, 100, "A100V", pp3_min_predictors = 1))
})
