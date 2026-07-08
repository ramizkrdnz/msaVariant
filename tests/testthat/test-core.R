test_that("build_msa_coord_map handles a simple alignment", {
  seqs <- c(ref = "AB-CD", other = "ABXCD")
  m <- build_msa_coord_map(seqs, ref_name = "ref")
  expect_equal(nrow(m), 5L)
  expect_equal(m$residue_pos, c(1L, 2L, NA_integer_, 3L, 4L))
  expect_equal(m$aa, c("A","B","-","C","D"))
})

test_that("map_variant_to_msa skips gaps in the reference", {
  seqs <- c(ref = "A-BC")
  expect_equal(map_variant_to_msa(c(1, 2, 3), seqs, "ref"), c(1L, 3L, 4L))
})

test_that("conservation_score is high for a perfectly conserved column", {
  seqs <- c(s1 = "AAAA", s2 = "AAAA", s3 = "AAAA")
  cs <- conservation_score(seqs, method = "shannon")
  expect_equal(cs$score, rep(1, 4))
})

test_that("conservation_score lowers for a fully variable column", {
  seqs <- c(s1 = "A", s2 = "C", s3 = "D", s4 = "E")
  cs <- conservation_score(seqs, method = "shannon")
  expect_true(cs$score[1] < 0.6)
})

test_that("conservation_score works with Jensen-Shannon method", {
  seqs <- c(s1 = "AAAA", s2 = "AAAA", s3 = "AAAA")
  cs <- conservation_score(seqs, method = "js")
  expect_true(all(!is.na(cs$score)))
})

test_that("geom_variant returns a list of layers", {
  seqs <- c(ref = "AAKAA", other = "AAKAA")
  v <- data.frame(pos = 3, label = "K3fs",
                  consequence = factor("frameshift",
                                       levels = c("frameshift","missense")))
  out <- geom_variant(v, msa = seqs, ref_name = "ref")
  expect_true(is.list(out))
  expect_true(length(out) >= 1L)
})

test_that("geom_variant errors when pos column missing", {
  seqs <- c(ref = "AAKAA")
  v <- data.frame(label = "K3fs")
  expect_error(geom_variant(v, msa = seqs, ref_name = "ref"),
               "must contain a `pos` column")
})

test_that("geom_domain renders a domain box", {
  seqs <- c(ref = "AAAAAAAAAA")
  d <- data.frame(start = 2, end = 7, name = "Test domain")
  out <- geom_domain(d, msa = seqs, ref_name = "ref")
  expect_true(is.list(out))
})

test_that("geom_track handles continuous data", {
  seqs <- c(ref = "AAAAA")
  d <- data.frame(pos = 1:5, score = c(0.1, 0.5, 0.9, 0.2, 0.7))
  out <- geom_track(d, msa = seqs, ref_name = "ref",
                    value = "score", type = "continuous",
                    name = "Demo")
  expect_true(is.list(out))
})

test_that("geom_track handles discrete data", {
  seqs <- c(ref = "AAAAA")
  d <- data.frame(pos = c(1, 3, 5),
                  sig = c("Pathogenic", "Benign", "VUS"))
  out <- geom_track(d, msa = seqs, ref_name = "ref",
                    value = "sig", type = "discrete",
                    name = "Demo")
  expect_true(is.list(out))
})

test_that("geom_track errors when value column missing", {
  seqs <- c(ref = "AAAAA")
  d <- data.frame(pos = 1:5, score = 1:5)
  expect_error(
    geom_track(d, msa = seqs, ref_name = "ref",
               value = "nonexistent", type = "continuous"),
    "not found in `data`"
  )
})

test_that("map_variant_to_msa warns on out-of-range positions", {
  seqs <- c(ref = "AAAA")
  expect_warning(
    map_variant_to_msa(c(1, 99), seqs, "ref"),
    "could not be mapped"
  )
})
