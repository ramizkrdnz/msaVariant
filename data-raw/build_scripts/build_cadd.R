## =================================================================
## build_cadd.R -- CADD v1.7 per-substitution scores
## =================================================================
##
## Reads the pre-VEP'd CADD inclAnno file. Extracts coding consequences
## and aggregates to per-substitution median PHRED.
##
## Storage: ~300 GB compressed. Streaming with `gzfile()` only;
## never read into memory.

build_cadd <- function(uniprot_ref, tmp_dir, european_mirror = FALSE) {

  base <- if (european_mirror)
    "https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh38"
  else
    "https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38"
  url <- file.path(base, "whole_genome_SNVs_inclAnno.tsv.gz")
  local <- file.path(tmp_dir, "cadd_inclAnno.tsv.gz")
  if (!file.exists(local)) {
    message("  Downloading CADD v1.7 inclAnno (~300 GB) ...")
    message("  This will take many hours.")
    utils::download.file(url, local, mode = "wb")
  } else {
    message("  Using cached CADD file")
  }

  cons_levels <- c("missense","synonymous","stop_gained",
                   "frameshift","inframe_deletion","inframe_insertion",
                   "splice_donor","splice_acceptor","other")

  # Map CADD consequence vocabulary to ours
  cons_map <- function(c) {
    case <- list(
      missense          = c("NON_SYNONYMOUS","MISSENSE"),
      synonymous        = c("SYNONYMOUS"),
      stop_gained       = c("STOP_GAINED"),
      frameshift        = c("FRAME_SHIFT"),
      inframe_deletion  = c("INFRAME_DELETION","DEL"),
      inframe_insertion = c("INFRAME_INSERTION","INS"),
      splice_donor      = c("CANONICAL_SPLICE","SPLICE_DONOR"),
      splice_acceptor   = c("SPLICE_ACCEPTOR")
    )
    for (k in names(case)) if (c %in% case[[k]]) return(k)
    "other"
  }

  message("  Streaming CADD inclAnno ...")
  con <- gzfile(local, "r")
  on.exit(close(con), add = TRUE)

  # First find which columns we need from the inclAnno header.
  # The header is ~120 columns; pick by name.
  header_lines <- readLines(con, n = 2L)
  # CADD has two header lines: ## comment + actual header
  hdr <- strsplit(header_lines[2], "\t", fixed = TRUE)[[1]]
  idx <- function(name) which(hdr == name)
  col <- list(
    gene  = idx("GeneName"),
    pp    = idx("ProteinPosition"),
    oAA   = idx("oAA"),
    nAA   = idx("nAA"),
    cons  = idx("Consequence"),
    raw   = idx("RawScore"),
    phred = idx("PHRED")
  )
  if (any(vapply(col, length, integer(1)) == 0L)) {
    stop("CADD header missing expected columns. ",
         "Inspect ", local, " and update the column names.")
  }

  buffers <- new.env(hash = TRUE, parent = emptyenv())
  n_rows <- 0L

  repeat {
    chunk <- readLines(con, n = 200000L)
    if (length(chunk) == 0L) break
    fields <- strsplit(chunk, "\t", fixed = TRUE)
    for (row in fields) {
      pp_s <- row[col$pp]
      if (is.na(pp_s) || pp_s == "NA" || !nzchar(pp_s)) next
      pos <- suppressWarnings(as.integer(pp_s))
      if (is.na(pos)) next
      gene <- row[col$gene]
      if (is.na(gene) || !nzchar(gene)) next
      oAA <- row[col$oAA]; nAA <- row[col$nAA]
      if (is.na(oAA) || nchar(oAA) != 1L) next
      aa_alt <- if (nchar(nAA) == 1L) nAA else tolower(nAA)
      consequence <- cons_map(row[col$cons])
      cadd_raw <- suppressWarnings(as.numeric(row[col$raw]))
      cadd_phred <- suppressWarnings(as.numeric(row[col$phred]))
      if (is.na(cadd_phred)) next
      aa_change <- paste0(oAA, pos, aa_alt)

      key <- paste0(gene, "::", aa_change)
      rec <- list(pos = pos, aa_ref = oAA, aa_alt = aa_alt,
                  aa_change = aa_change, consequence = consequence,
                  cadd_raw = cadd_raw, cadd_phred = cadd_phred,
                  gene = gene)

      if (is.null(buffers[[key]])) buffers[[key]] <- list()
      buffers[[key]][[length(buffers[[key]]) + 1L]] <- rec
      n_rows <- n_rows + 1L
    }
  }
  message(sprintf("  Parsed %d CADD rows aggregating to %d unique substitutions",
                  n_rows, length(ls(buffers))))

  # Aggregate by aa_change taking median PHRED
  per_gene <- list()
  for (key in ls(buffers)) {
    rows <- buffers[[key]]
    first <- rows[[1]]
    gene <- first$gene
    rec <- data.frame(
      pos         = first$pos,
      aa_ref      = first$aa_ref,
      aa_alt      = first$aa_alt,
      aa_change   = first$aa_change,
      consequence = first$consequence,
      cadd_raw    = stats::median(vapply(rows, `[[`, numeric(1), "cadd_raw")),
      cadd_phred  = stats::median(vapply(rows, `[[`, numeric(1), "cadd_phred")),
      stringsAsFactors = FALSE
    )
    if (is.null(per_gene[[gene]])) per_gene[[gene]] <- list()
    per_gene[[gene]][[length(per_gene[[gene]]) + 1L]] <- rec
  }

  out <- lapply(per_gene, function(rows) {
    df <- do.call(rbind, rows)
    df$pos <- as.integer(df$pos)
    df$consequence <- factor(df$consequence, levels = cons_levels)
    rownames(df) <- NULL
    df
  })
  out
}
