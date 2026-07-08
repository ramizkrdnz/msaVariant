## =================================================================
## build_revel.R -- REVEL per-substitution scores
## =================================================================
##
## REVEL is keyed by genomic position. The published file
## (revel-v1.3_all_chromosomes.zip from Mt. Sinai) contains an
## `aapos` column when present, otherwise we need VEP annotation
## first. Verify the file's column layout at build time.

build_revel <- function(uniprot_ref, tmp_dir) {

  url <- "https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip"
  local <- file.path(tmp_dir, "revel.zip")
  if (!file.exists(local)) {
    message("  Downloading REVEL (~700 MB) ...")
    utils::download.file(url, local, mode = "wb")
  }
  unzdir <- file.path(tmp_dir, "revel_unz")
  if (!dir.exists(unzdir)) {
    message("  Unzipping REVEL ...")
    utils::unzip(local, exdir = unzdir)
  }
  # The zip contains one CSV (sometimes chromosome-split); find it
  csvs <- list.files(unzdir, pattern = "\\.csv$",
                      recursive = TRUE, full.names = TRUE)
  if (length(csvs) == 0L) stop("No CSV files found inside REVEL zip.")

  ens2g <- setNames(uniprot_ref$gene,
                    uniprot_ref$ensembl_transcript_id)
  ens2g <- ens2g[!is.na(names(ens2g))]

  message("  Streaming REVEL CSV(s) ...")
  out <- list()
  for (csv in csvs) {
    df <- readr::read_csv(csv, show_col_types = FALSE, progress = FALSE)
    # Tolerant column lookup
    pick <- function(p) {
      m <- grep(p, names(df), ignore.case = TRUE, value = TRUE)[1]
      if (is.na(m)) return(NULL)
      df[[m]]
    }
    aaref <- pick("^aaref$"); aaalt <- pick("^aaalt$")
    aapos <- pick("^aapos$|^protein.position$|^aa.pos$")
    revel <- pick("^REVEL$|^revel$|^revel_score$")
    tx    <- pick("ensembl.*transcript|^transcript$")
    if (any(sapply(list(aaref, aaalt, aapos, revel, tx), is.null))) {
      warning("REVEL file ", csv, " is missing required columns; skipping.")
      next
    }

    # Multi-transcript rows: REVEL has semicolon-separated lists
    # in tx and aapos. Split and flatten.
    # For simplicity, take first transcript per row.
    first <- function(x) sub(";.*", "", x)
    sub_df <- data.frame(
      aa_ref      = aaref,
      aa_alt      = aaalt,
      pos         = suppressWarnings(as.integer(first(as.character(aapos)))),
      revel_score = as.numeric(revel),
      tx          = first(as.character(tx)),
      stringsAsFactors = FALSE
    )
    sub_df <- sub_df[!is.na(sub_df$pos) & !is.na(sub_df$revel_score), ]
    sub_df$gene <- unname(ens2g[sub_df$tx])
    sub_df <- sub_df[!is.na(sub_df$gene), ]
    sub_df$aa_change <- paste0(sub_df$aa_ref, sub_df$pos, sub_df$aa_alt)

    # Slice and aggregate per gene + aa_change (max REVEL across
    # codon-redundant SNVs)
    for (g in unique(sub_df$gene)) {
      gd <- sub_df[sub_df$gene == g, ]
      agg <- stats::aggregate(revel_score ~ aa_ref + aa_alt + pos + aa_change,
                              data = gd, FUN = max)
      agg$pos <- as.integer(agg$pos)
      agg <- agg[, c("pos", "aa_ref", "aa_alt", "aa_change", "revel_score")]
      rownames(agg) <- NULL
      if (g %in% names(out)) {
        out[[g]] <- unique(rbind(out[[g]], agg))
      } else {
        out[[g]] <- agg
      }
    }
  }
  message(sprintf("  Wrote %d genes", length(out)))
  out
}
