## =================================================================
## build_alphamissense.R -- AlphaMissense_aa_substitutions.tsv.gz
## =================================================================

build_alphamissense <- function(uniprot_ref, tmp_dir) {

  url <- "https://zenodo.org/record/8208688/files/AlphaMissense_aa_substitutions.tsv.gz"
  local <- file.path(tmp_dir, "AlphaMissense_aa_substitutions.tsv.gz")
  if (!file.exists(local)) {
    message("  Downloading AlphaMissense (1.2 GB) ...")
    utils::download.file(url, local, mode = "wb")
    # Verify checksum
    obs_md5 <- tools::md5sum(local)
    exp_md5 <- "b9ccb339e0de6cb0a8d1973ad2026576"
    if (unname(obs_md5) != exp_md5) {
      stop("MD5 mismatch on AlphaMissense download. ",
           "Got ", obs_md5, " expected ", exp_md5,
           ". Delete the file and retry.")
    }
  } else {
    message("  Using cached AlphaMissense file")
  }

  u2g <- setNames(uniprot_ref$gene, uniprot_ref$uniprot_id)
  am_levels <- c("likely_benign","ambiguous","likely_pathogenic")

  message("  Streaming AlphaMissense TSV ...")
  con <- gzfile(local, "r")
  on.exit(close(con), add = TRUE)

  # Skip header comment lines (start with #), then header row
  repeat {
    ln <- readLines(con, n = 1L)
    if (length(ln) == 0L) break
    if (!startsWith(ln, "#")) break  # this is the header row
  }
  # The non-# header line was just read; ln is the header row.

  buffers <- new.env(hash = TRUE, parent = emptyenv())
  n_rows <- 0L

  repeat {
    chunk <- readLines(con, n = 500000L)
    if (length(chunk) == 0L) break
    fields <- strsplit(chunk, "\t", fixed = TRUE)
    for (row in fields) {
      if (length(row) < 4L) next
      acc <- sub("-\\d+$", "", row[1])
      gene <- u2g[[acc]]
      if (is.null(gene) || is.na(gene)) next
      variant <- row[2]  # e.g. "K518R"
      score   <- suppressWarnings(as.numeric(row[3]))
      cls     <- row[4]
      pos <- suppressWarnings(as.integer(sub("^[A-Z](\\d+)[A-Z]$", "\\1", variant)))
      if (is.na(pos) || is.na(score)) next
      rec <- c(pos = pos,
               aa_ref = substr(variant, 1, 1),
               aa_alt = substr(variant, nchar(variant), nchar(variant)),
               aa_change = variant,
               am_score = score,
               am_class = cls)
      key <- gene
      if (is.null(buffers[[key]])) buffers[[key]] <- list()
      buffers[[key]][[length(buffers[[key]]) + 1L]] <- rec
      n_rows <- n_rows + 1L
    }
  }
  message(sprintf("  Parsed %d AlphaMissense substitutions across %d genes",
                  n_rows, length(ls(buffers))))

  out <- list()
  for (gene in ls(buffers)) {
    rows <- buffers[[gene]]
    df <- data.frame(
      pos       = as.integer(vapply(rows, `[[`, character(1), "pos")),
      aa_ref    = vapply(rows, `[[`, character(1), "aa_ref"),
      aa_alt    = vapply(rows, `[[`, character(1), "aa_alt"),
      aa_change = vapply(rows, `[[`, character(1), "aa_change"),
      am_score  = as.numeric(vapply(rows, `[[`, character(1), "am_score")),
      am_class  = factor(vapply(rows, `[[`, character(1), "am_class"),
                         levels = am_levels),
      stringsAsFactors = FALSE
    )
    rownames(df) <- NULL
    out[[gene]] <- df
  }
  out
}
