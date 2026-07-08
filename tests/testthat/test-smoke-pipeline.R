## End-to-end pipeline smoke test.
##
## Converted from the former standalone tests/smoke_test.R (which had
## hard-coded /home/claude/... paths and failed R CMD check). It drives
## the whole PATL1 ortholog pipeline -- coord map, geoms, conservation,
## generic track, and a full ggplot assembly saved to a temp dir -- so
## the integration path is exercised, not just the individual units.

fa <- system.file("extdata", "patl1_orthologs.fasta",
                  package = "msaVariant")

test_that("PATL1 ortholog fixture is available", {
  expect_true(nzchar(fa) && file.exists(fa))
})

test_that("build_msa_coord_map maps the K518 column", {
  skip_if(!nzchar(fa) || !file.exists(fa))
  m <- build_msa_coord_map(fa, ref_name = "PATL1_HUMAN")
  # column 50 corresponds to K518 in PATL1 numbering
  expect_identical(m$aa[50], "K")
})

test_that("conservation_score finds K518 invariant", {
  skip_if(!nzchar(fa) || !file.exists(fa))
  cons <- conservation_score(fa, method = "shannon")
  k518_score <- cons$score[cons$msa_col == 50]
  expect_gt(k518_score, 0.9)
})

test_that("geom_variant / geom_domain / geom_track return ggplot layers", {
  skip_if(!nzchar(fa) || !file.exists(fa))
  m <- build_msa_coord_map(fa, ref_name = "PATL1_HUMAN")
  res_map <- m[!is.na(m$residue_pos), ]

  variants <- data.frame(
    pos         = 518 - 490 + 1,                # = 29
    pos_end     = max(res_map$residue_pos),
    label       = "K518fs",
    consequence = factor("frameshift",
                         levels = c("frameshift", "nonsense", "missense"))
  )
  v_layers <- geom_variant(variants, msa = fa, ref_name = "PATL1_HUMAN")
  expect_true(is.list(v_layers))
  expect_gt(length(v_layers), 0L)

  domains <- data.frame(
    start = 1, end = max(res_map$residue_pos),
    name  = "PAT1 mid-domain (LSm1-7 binding region)"
  )
  dom_layers <- geom_domain(domains, msa = fa, ref_name = "PATL1_HUMAN")
  expect_true(is.list(dom_layers))

  cons <- conservation_score(fa, method = "shannon")
  cons_df <- data.frame(
    pos          = m$residue_pos[!is.na(m$residue_pos)],
    conservation = cons$score[!is.na(m$residue_pos)]
  )
  cons_df <- cons_df[!is.na(cons_df$conservation), ]
  trk_cons <- geom_track(cons_df, msa = fa, ref_name = "PATL1_HUMAN",
                         value = "conservation", type = "continuous",
                         name  = "Conservation",
                         value_range = c(0, 1),
                         palette = c("#FFFFFF", "#FFE082", "#FB8C00", "#B71C1C"),
                         y_offset = -2, track_height = 1.0)
  expect_false(is.null(trk_cons))
})

test_that("full MSA + variant + domain + track plot assembles and saves", {
  skip_if(!nzchar(fa) || !file.exists(fa))
  m <- build_msa_coord_map(fa, ref_name = "PATL1_HUMAN")
  res_map <- m[!is.na(m$residue_pos), ]

  variants <- data.frame(
    pos         = 518 - 490 + 1,
    pos_end     = max(res_map$residue_pos),
    label       = "K518fs",
    consequence = factor("frameshift",
                         levels = c("frameshift", "nonsense", "missense"))
  )
  v_layers   <- geom_variant(variants, msa = fa, ref_name = "PATL1_HUMAN")
  domains    <- data.frame(start = 1, end = max(res_map$residue_pos),
                           name = "PAT1 mid-domain")
  dom_layers <- geom_domain(domains, msa = fa, ref_name = "PATL1_HUMAN")

  cons <- conservation_score(fa, method = "shannon")
  cons_df <- data.frame(pos = m$residue_pos[!is.na(m$residue_pos)],
                        conservation = cons$score[!is.na(m$residue_pos)])
  cons_df <- cons_df[!is.na(cons_df$conservation), ]
  trk_cons <- geom_track(cons_df, msa = fa, ref_name = "PATL1_HUMAN",
                         value = "conservation", type = "continuous",
                         name = "Conservation", value_range = c(0, 1),
                         palette = c("#FFFFFF", "#FFE082", "#FB8C00", "#B71C1C"),
                         y_offset = -2, track_height = 1.0)

  seqs <- {
    ls <- readLines(fa); is_h <- startsWith(ls, ">")
    hdr <- sub("^>", "", ls[is_h])
    stats::setNames(ls[!is_h], hdr)
  }
  mat_df <- do.call(rbind, lapply(seq_along(seqs), function(i) {
    cs <- strsplit(seqs[i], "", fixed = TRUE)[[1]]
    data.frame(seq_name = names(seqs)[i], col = seq_along(cs),
               aa = cs, row = i, stringsAsFactors = FALSE)
  }))

  p <- ggplot2::ggplot(mat_df) +
    ggplot2::geom_tile(ggplot2::aes(x = col, y = seq_name, fill = aa),
                       colour = "white", linewidth = 0.2) +
    v_layers + dom_layers + trk_cons +
    ggplot2::labs(x = "PATL1 residue position", y = NULL)

  expect_s3_class(p, "ggplot")

  # Outputs go to a self-cleaning temp dir, never a hard-coded path.
  out_dir <- withr::local_tempdir()
  png_path <- file.path(out_dir, "msaVariant_demo.png")
  ggplot2::ggsave(png_path, p, width = 16, height = 5.5, dpi = 90)
  expect_true(file.exists(png_path))
  expect_gt(file.info(png_path)$size, 0L)
})
