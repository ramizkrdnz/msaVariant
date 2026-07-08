## =================================================================
## build_uniprot_ref.R
## =================================================================
##
## Builds the master UniProt -> HGNC -> Ensembl reference for all
## reviewed (SwissProt) human protein entries (~20,000 rows).
##
## Returns a data.frame with columns:
##   uniprot_id, gene, gene_synonyms, protein_length,
##   ensembl_gene_id, ensembl_transcript_id, ensembl_protein_id,
##   refseq_protein, mane_select, sequence_version
##
## Cached to disk on first run so subsequent rebuilds skip the
## ~5 MB download.

build_uniprot_ref <- function(tmp_dir) {
  cache <- file.path(tmp_dir, "uniprot_ref.rds")
  if (file.exists(cache)) {
    return(readRDS(cache))
  }

  url <- paste0(
    "https://rest.uniprot.org/uniprotkb/stream?",
    "query=organism_id:9606+AND+reviewed:true&",
    "format=tsv&",
    "fields=accession,gene_primary,gene_names,length,",
    "sequence_version,xref_ensembl,xref_refseq,xref_mane-select"
  )
  raw <- file.path(tmp_dir, "uniprot_human_ref.tsv")
  if (!file.exists(raw)) {
    message("  Downloading UniProt human reference (~5 MB) ...")
    utils::download.file(url, raw, mode = "wb")
  }

  df <- readr::read_tsv(raw, show_col_types = FALSE)
  # UniProt's TSV column names are inconsistent; normalize them
  names(df) <- gsub("\\.\\.\\.", "_", make.names(names(df)))

  # Tolerant column lookup: pick the first column whose name
  # matches each candidate pattern.
  pick <- function(patterns) {
    for (p in patterns) {
      m <- grep(p, names(df), ignore.case = TRUE, value = TRUE)
      if (length(m)) return(df[[m[1]]])
    }
    rep(NA_character_, nrow(df))
  }

  out <- data.frame(
    uniprot_id            = pick("^Entry$|accession"),
    gene                  = pick("Gene.*Primary|gene_primary"),
    gene_synonyms         = pick("Gene.*Names|gene_names"),
    protein_length        = as.integer(pick("^Length$|length")),
    ensembl_gene_id       = .extract_ensembl_id(pick("Ensembl"),
                                                  prefix = "ENSG"),
    ensembl_transcript_id = .extract_ensembl_id(pick("Ensembl"),
                                                  prefix = "ENST"),
    ensembl_protein_id    = .extract_ensembl_id(pick("Ensembl"),
                                                  prefix = "ENSP"),
    refseq_protein        = .extract_first(pick("RefSeq"), "^NP_"),
    mane_select           = .extract_first(pick("MANE"), "^NM_"),
    sequence_version      = as.integer(pick("Sequence.*Version|sequence_version")),
    stringsAsFactors      = FALSE
  )

  # Take the first primary gene name (some entries have multiple)
  out$gene <- sub("\\s.*", "", out$gene)
  # Strip isoform suffixes from uniprot IDs (we want canonical only)
  out$uniprot_id <- sub("-\\d+$", "", out$uniprot_id)
  # Drop entries without an HGNC symbol or with bad length
  out <- out[!is.na(out$gene) & nzchar(out$gene) &
              !is.na(out$protein_length) & out$protein_length > 0L, ]

  saveRDS(out, cache)
  out
}

# Helper: extract first instance of a prefix from a "key1;key2;..." string
.extract_ensembl_id <- function(x, prefix) {
  out <- character(length(x))
  for (i in seq_along(x)) {
    s <- x[i]
    if (is.na(s) || !nzchar(s)) { out[i] <- NA_character_; next }
    ids <- strsplit(s, "[;,\\s]+")[[1]]
    hit <- grep(paste0("^", prefix), ids, value = TRUE)
    out[i] <- if (length(hit)) hit[1] else NA_character_
  }
  out
}

.extract_first <- function(x, pattern) {
  out <- character(length(x))
  for (i in seq_along(x)) {
    s <- x[i]
    if (is.na(s) || !nzchar(s)) { out[i] <- NA_character_; next }
    ids <- strsplit(s, "[;,\\s]+")[[1]]
    hit <- grep(pattern, ids, value = TRUE)
    out[i] <- if (length(hit)) hit[1] else NA_character_
  }
  out
}
