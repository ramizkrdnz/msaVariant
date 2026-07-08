## =================================================================
## build_domains.R -- InterPro protein2ipr.dat.gz slicer
## =================================================================
##
## Reads protein2ipr.dat.gz, joins to UniProt -> HGNC, produces a
## named list (gene -> domains data.frame) conforming to the
## DATA_FORMAT_SPEC.md `domains` schema.

build_domains <- function(uniprot_ref, tmp_dir) {

  url <- "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/protein2ipr.dat.gz"
  local <- file.path(tmp_dir, "protein2ipr.dat.gz")
  if (!file.exists(local)) {
    message("  Downloading InterPro protein2ipr.dat.gz (~5 GB) ...")
    utils::download.file(url, local, mode = "wb")
  } else {
    message("  Using cached protein2ipr.dat.gz")
  }

  # Build UniProt -> gene lookup vector
  u2g <- setNames(uniprot_ref$gene, uniprot_ref$uniprot_id)

  # Member-database -> canonical source label map
  # Pfam: PF00000, InterPro: IPR000000, SMART: SM00000, etc.
  source_for <- function(sig) {
    if (grepl("^PF", sig)) "Pfam"
    else if (grepl("^IPR", sig)) "InterPro"
    else if (grepl("^SM", sig)) "SMART"
    else if (grepl("^PS", sig)) "PROSITE"
    else if (grepl("^PR", sig)) "PRINTS"
    else if (grepl("^PTHR", sig)) "PANTHER"
    else NA_character_
  }
  source_levels <- c("Pfam", "InterPro", "SMART", "PROSITE",
                     "PRINTS", "PANTHER")

  message("  Streaming protein2ipr.dat.gz ...")
  con <- gzfile(local, "r")
  on.exit(close(con), add = TRUE)

  buffers <- list()
  n_rows <- 0L

  repeat {
    chunk <- readLines(con, n = 200000L)
    if (length(chunk) == 0L) break
    fields <- strsplit(chunk, "\t", fixed = TRUE)
    for (row in fields) {
      if (length(row) < 6L) next
      acc  <- sub("-\\d+$", "", row[1])  # strip isoform suffix
      gene <- u2g[[acc]]
      if (is.null(gene) || is.na(gene)) next
      sig <- row[4]
      src <- source_for(sig)
      if (is.na(src)) next
      rec <- list(
        start     = as.integer(row[5]),
        end       = as.integer(row[6]),
        name      = row[3],
        accession = row[2],
        source    = src
      )
      if (is.null(buffers[[gene]])) buffers[[gene]] <- list()
      buffers[[gene]][[length(buffers[[gene]]) + 1L]] <- rec
      n_rows <- n_rows + 1L
    }
  }
  message(sprintf("  Parsed %d domain records across %d genes",
                  n_rows, length(buffers)))

  # Convert buffers to per-gene data.frames
  out <- lapply(buffers, function(records) {
    df <- do.call(rbind.data.frame, c(records, stringsAsFactors = FALSE))
    df <- df[order(df$start, df$end), , drop = FALSE]
    df$source <- factor(df$source, levels = source_levels)
    rownames(df) <- NULL
    df
  })
  out
}
