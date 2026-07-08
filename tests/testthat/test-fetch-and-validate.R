# Helper: build a valid synthetic gene bundle for validator tests
.make_valid_bundle <- function(gene = "TESTGENE",
                                protein_length = 100L) {
  list(
    meta = data.frame(
      gene                   = gene,
      uniprot_id             = "P00000",
      protein_length         = protein_length,
      ensembl_gene_id        = "ENSG00000000000",
      ensembl_transcript_id  = "ENST00000000000",
      build_date             = Sys.Date(),
      source_versions        = "ClinVar=2025-01;gnomAD=v4.1",
      stringsAsFactors       = FALSE
    ),
    domains = data.frame(
      start     = 10L,
      end       = 50L,
      name      = "Test domain",
      accession = "PF00000",
      source    = factor("Pfam",
                         levels = c("Pfam","InterPro","SMART",
                                    "PROSITE","PRINTS","PANTHER")),
      stringsAsFactors = FALSE
    ),
    clinvar = data.frame(
      pos          = 25L,
      aa_ref       = "K",
      aa_alt       = "R",
      aa_change    = "K25R",
      significance = factor("Pathogenic",
                            levels = c("Pathogenic","Likely_pathogenic","VUS",
                                       "Likely_benign","Benign","Conflicting","Other")),
      review_status = factor("3_star",
                             levels = c("4_star","3_star","2_star","1_star","0_star")),
      clinvar_id   = "VCV000001",
      stringsAsFactors = FALSE
    ),
    gnomad = data.frame(
      pos          = 30L,
      aa_ref       = "A",
      aa_alt       = "V",
      aa_change    = "A30V",
      consequence  = factor("missense",
                            levels = c("missense","synonymous","stop_gained",
                                       "frameshift","inframe_deletion","inframe_insertion",
                                       "splice_donor","splice_acceptor","other")),
      af_joint     = 1e-5,
      ac_joint     = 1L,
      an_joint     = 100000L,
      filter       = "PASS",
      stringsAsFactors = FALSE
    ),
    alphamissense = data.frame(
      pos       = 25L,
      aa_ref    = "K",
      aa_alt    = "R",
      aa_change = "K25R",
      am_score  = 0.85,
      am_class  = factor("likely_pathogenic",
                         levels = c("likely_benign","ambiguous","likely_pathogenic")),
      stringsAsFactors = FALSE
    ),
    revel = data.frame(
      pos         = 25L,
      aa_ref      = "K",
      aa_alt      = "R",
      aa_change   = "K25R",
      revel_score = 0.72,
      stringsAsFactors = FALSE
    ),
    cadd = data.frame(
      pos         = 25L,
      aa_ref      = "K",
      aa_alt      = "R",
      aa_change   = "K25R",
      consequence = factor("missense",
                           levels = c("missense","synonymous","stop_gained",
                                      "frameshift","inframe_deletion","inframe_insertion",
                                      "splice_donor","splice_acceptor","other")),
      cadd_raw    = 4.2,
      cadd_phred  = 28.5,
      stringsAsFactors = FALSE
    )
  )
}

test_that("validate_gene_data accepts a valid bundle", {
  b <- .make_valid_bundle()
  v <- validate_gene_data(b, strict = FALSE)
  expect_true(v$valid)
  expect_length(v$issues, 0L)
})

test_that("validate_gene_data rejects missing top-level element", {
  b <- .make_valid_bundle()
  b$clinvar <- NULL
  v <- validate_gene_data(b, strict = FALSE)
  expect_false(v$valid)
  expect_match(v$issues[1], "Missing top-level element")
})

test_that("validate_gene_data rejects wrong column type", {
  b <- .make_valid_bundle()
  b$gnomad$af_joint <- as.character(b$gnomad$af_joint)
  v <- validate_gene_data(b, strict = FALSE)
  expect_false(v$valid)
  expect_match(paste(v$issues, collapse = " | "), "af_joint")
})

test_that("validate_gene_data rejects wrong factor levels", {
  b <- .make_valid_bundle()
  b$clinvar$significance <- factor("Pathogenic")  # one-level factor
  v <- validate_gene_data(b, strict = FALSE)
  expect_false(v$valid)
  expect_match(paste(v$issues, collapse = " | "), "levels mismatch")
})

test_that("validate_gene_data rejects NA in required columns", {
  b <- .make_valid_bundle()
  b$clinvar$aa_change[1] <- NA_character_
  v <- validate_gene_data(b, strict = FALSE)
  expect_false(v$valid)
  expect_match(paste(v$issues, collapse = " | "), "NA in required")
})

test_that("validate_gene_data flags out-of-range positions", {
  b <- .make_valid_bundle(protein_length = 50L)
  b$clinvar$pos <- 9999L
  v <- validate_gene_data(b, strict = FALSE)
  expect_false(v$valid)
  expect_match(paste(v$issues, collapse = " | "), "outside")
})

test_that("validate_gene_data accepts empty tables", {
  b <- .make_valid_bundle()
  b$clinvar <- b$clinvar[0, ]   # empty but with right columns
  v <- validate_gene_data(b, strict = FALSE)
  expect_true(v$valid)
})

test_that("validate_gene_data strict mode aborts on failure", {
  b <- .make_valid_bundle()
  b$meta <- NULL
  expect_error(validate_gene_data(b, strict = TRUE), "failed validation")
})

test_that("fetch_gene_data aborts cleanly when Zenodo URL is unset", {
  expect_error(fetch_gene_data("PATL1"),
               "has not been configured")
})

test_that("fetch_gene_data rejects suspicious gene symbols", {
  expect_error(fetch_gene_data("../../etc/passwd"),
               "Suspicious gene symbol")
})

test_that("cache_summary works on empty cache", {
  # Use a temp cache dir to avoid clobbering real cache
  old <- Sys.getenv("MSAVARIANT_CACHE")
  on.exit(Sys.setenv(MSAVARIANT_CACHE = old))
  Sys.setenv(MSAVARIANT_CACHE = tempfile())
  s <- cache_summary()
  expect_s3_class(s, "data.frame")
  expect_equal(nrow(s), 0L)
})

test_that("geom_gnomad accepts user data without network", {
  seqs <- c(ref = "AAAAA")
  d <- data.frame(pos = c(1, 3),
                  af_joint = c(1e-4, 5e-5),
                  stringsAsFactors = FALSE)
  out <- geom_gnomad(data = d, msa = seqs, ref_name = "ref")
  expect_true(is.list(out))
})

test_that("geom_clinvar accepts user data without network", {
  seqs <- c(ref = "AAAAA")
  d <- data.frame(
    pos = c(1, 3),
    significance = factor(c("Pathogenic", "Benign"),
                          levels = c("Pathogenic","Likely_pathogenic","VUS",
                                     "Likely_benign","Benign","Conflicting","Other")),
    stringsAsFactors = FALSE
  )
  out <- geom_clinvar(data = d, msa = seqs, ref_name = "ref")
  expect_true(is.list(out))
})

test_that("geom_alphamissense aggregates per-substitution to per-residue", {
  seqs <- c(ref = "AAAAA")
  # 3 alts at residue 2, plus one at residue 4
  d <- data.frame(
    pos = c(2, 2, 2, 4),
    am_score = c(0.1, 0.9, 0.5, 0.7),
    stringsAsFactors = FALSE
  )
  out <- geom_alphamissense(data = d, msa = seqs, ref_name = "ref")
  expect_true(is.list(out))
})

test_that("geom_* error when both data and gene supplied", {
  seqs <- c(ref = "AAAA")
  d <- data.frame(pos = 1, af_joint = 0.001)
  expect_error(geom_gnomad(data = d, gene = "TP53",
                            msa = seqs, ref_name = "ref"),
               "pass either")
})
