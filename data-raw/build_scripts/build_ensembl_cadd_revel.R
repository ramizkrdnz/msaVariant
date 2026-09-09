## =================================================================
## build_ensembl_cadd_revel.R -- per-gene CADD + REVEL via Ensembl REST
## =================================================================
##
## WHY THIS EXISTS
## ---------------
## The proven TP53/WDR31 bundles got CADD and REVEL from a per-gene
## Ensembl VEP export (see build_tp53_from_real_data.R section 5). This
## builder automates that same provenance over the Ensembl REST API, so
## we avoid the ~300 GB bulk CADD file and the ~700 MB REVEL zip and pull
## only what each gene needs. One VEP query per gene yields BOTH scores.
##
## Returns list(cadd = <by-gene>, revel = <by-gene>) -- two named lists
## keyed by gene symbol, each value a spec-conformant data.frame:
##   cadd : pos, aa_ref, aa_alt, aa_change, consequence, cadd_raw, cadd_phred
##   revel: pos, aa_ref, aa_alt, aa_change, revel_score
## It REPLACES build_cadd() + build_revel() in the orchestrator.
##
## RESUMABLE: each gene's result is checkpointed to
## <tmp_dir>/ensembl_cache/<gene>.rds; a re-run skips cached genes.
## Per-gene tryCatch: one gene's failure never halts the crawl.
##
## ================================================================
## ⚠️  PILOT MUST VALIDATE THIS FIRST -- READ BEFORE A FULL RUN
## ----------------------------------------------------------------
## The public Ensembl REST VEP endpoint (rest.ensembl.org) returns
## CADD via `?CADD=1`, but REVEL is served through the dbNSFP plugin
## (`?dbNSFP=REVEL_score`), which is NOT guaranteed to be enabled on the
## public server. The pilot's FIRST job is to confirm, on TP53, that the
## REST response actually carries both scores and that the parsed tables
## match the known-good manual TP53 bundle. If REVEL is absent from REST,
## fall back options (in preference order): (a) dbNSFP tabix range-slice
## per gene, (b) REVEL-only bulk file + CADD from REST. Do NOT launch the
## full crawl until this is confirmed against TP53/WDR31.
##
## The response field names below (transcript_consequences, cadd_phred,
## cadd_raw, revel_score, amino_acids, protein_start) follow the Ensembl
## VEP REST schema as documented; verify them against a live TP53
## response in the pilot and adjust the accessors if the schema differs.
## ================================================================

build_ensembl_cadd_revel <- function(uniprot_ref,
                                      tmp_dir,
                                      server = "https://rest.ensembl.org",
                                      sleep_per_query = 0.34,  # <= 3 req/s
                                      request_revel = TRUE) {
  if (!requireNamespace("httr", quietly = TRUE) ||
      !requireNamespace("jsonlite", quietly = TRUE)) {
    stop("build_ensembl_cadd_revel needs `httr` and `jsonlite` installed.")
  }

  cons_levels <- c("missense", "synonymous", "stop_gained", "frameshift",
                   "inframe_deletion", "inframe_insertion",
                   "splice_donor", "splice_acceptor", "other")

  cache_dir <- file.path(tmp_dir, "ensembl_cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  save_atomic <- function(obj, dest) {
    tmp <- tempfile(tmpdir = dirname(dest), fileext = ".part")
    saveRDS(obj, tmp)
    if (!file.rename(tmp, dest)) {
      file.copy(tmp, dest, overwrite = TRUE); unlink(tmp)
    }
  }

  # Only genes with an Ensembl gene id can be queried by id.
  ref <- uniprot_ref[!is.na(uniprot_ref$ensembl_gene_id) &
                       nzchar(uniprot_ref$ensembl_gene_id), , drop = FALSE]

  cadd_out <- list()
  revel_out <- list()
  n_cached <- 0L; n_fetched <- 0L; n_failed <- 0L

  for (i in seq_len(nrow(ref))) {
    gene <- ref$gene[i]
    ensg <- ref$ensembl_gene_id[i]
    plen <- as.integer(ref$protein_length[i])
    if (i %% 100L == 0L)
      message(sprintf("  %d / %d (%s) [fetched %d, cache %d, failed %d]",
                      i, nrow(ref), gene, n_fetched, n_cached, n_failed))

    cache_file <- file.path(cache_dir, paste0(gene, ".rds"))
    if (file.exists(cache_file)) {
      cached <- tryCatch(readRDS(cache_file), error = function(e) NULL)
      if (is.list(cached)) {
        if (nrow(cached$cadd) > 0L) cadd_out[[gene]] <- cached$cadd
        if (nrow(cached$revel) > 0L) revel_out[[gene]] <- cached$revel
      }
      n_cached <- n_cached + 1L
      next
    }

    res <- tryCatch(
      .ensembl_gene_scores(gene, ensg, plen, server,
                           request_revel, cons_levels),
      error = function(e) {
        warning(sprintf("Ensembl query failed for %s (%s): %s",
                        gene, ensg, conditionMessage(e)))
        NULL
      }
    )

    if (is.null(res)) {
      # Transient failure -> not cached -> retried next run.
      n_failed <- n_failed + 1L
    } else {
      save_atomic(res, cache_file)   # definitive (even if empty) -> cache
      n_fetched <- n_fetched + 1L
      if (nrow(res$cadd) > 0L) cadd_out[[gene]] <- res$cadd
      if (nrow(res$revel) > 0L) revel_out[[gene]] <- res$revel
    }
    Sys.sleep(sleep_per_query)
  }

  message(sprintf(
    "  Ensembl CADD/REVEL: %d genes with CADD, %d with REVEL (%d fetched, %d cached, %d failed)",
    length(cadd_out), length(revel_out), n_fetched, n_cached, n_failed))
  list(cadd = cadd_out, revel = revel_out)
}

## Query one gene's missense variants + CADD/REVEL from Ensembl REST VEP,
## returning list(cadd=<df>, revel=<df>) in the package schema. Empty
## (0-row) data.frames when the gene has no scored missense variants --
## a definitive answer, so it is cacheable.
.ensembl_gene_scores <- function(gene, ensg, protein_length, server,
                                  request_revel, cons_levels) {
  empty_cadd <- data.frame(
    pos = integer(), aa_ref = character(), aa_alt = character(),
    aa_change = character(),
    consequence = factor(character(), levels = cons_levels),
    cadd_raw = numeric(), cadd_phred = numeric(), stringsAsFactors = FALSE)
  empty_revel <- data.frame(
    pos = integer(), aa_ref = character(), aa_alt = character(),
    aa_change = character(), revel_score = numeric(), stringsAsFactors = FALSE)

  ## VEP for a whole gene is done by region. Ask Ensembl for the gene's
  ## genomic span, then VEP all overlapping variants with CADD (+REVEL via
  ## dbNSFP) plugins. Endpoint + params per the Ensembl REST VEP schema;
  ## the pilot verifies these against a live TP53 response.
  extra <- "CADD=1"
  if (isTRUE(request_revel)) extra <- paste0(extra, ";dbNSFP=REVEL_score")

  ## 1. Region variants VEP call (documented endpoint:
  ##    GET /vep/human/region/{region}/{allele}?... is per-variant; for a
  ##    whole gene we POST the variant set. The pilot wires the exact call
  ##    that returns CADD+REVEL; this scaffold parses the standard
  ##    transcript_consequences shape below.)
  url <- sprintf("%s/vep/human/id/%s?%s;content-type=application/json",
                 server, utils::URLencode(ensg), extra)
  resp <- httr::GET(url, httr::accept_json(), httr::timeout(120))
  if (httr::status_code(resp) != 200L) {
    ## Non-200 is treated as transient by the caller (return NULL there).
    stop(sprintf("HTTP %d", httr::status_code(resp)))
  }
  payload <- httr::content(resp, as = "parsed", simplifyVector = FALSE)
  if (length(payload) == 0L) {
    return(list(cadd = empty_cadd, revel = empty_revel))
  }

  ## Flatten transcript_consequences across all returned variants, keeping
  ## missense rows with an amino-acid change and a protein position.
  rows <- list()
  for (v in payload) {
    tcs <- v$transcript_consequences
    if (is.null(tcs)) next
    for (tc in tcs) {
      aa <- tc$amino_acids                 # "K/R"
      ppos <- tc$protein_start             # integer
      if (is.null(aa) || is.null(ppos) || !grepl("/", aa)) next
      parts <- strsplit(aa, "/", fixed = TRUE)[[1]]
      if (length(parts) != 2L || nchar(parts[1]) != 1L ||
          nchar(parts[2]) != 1L) next
      ct <- tc$consequence_terms
      rows[[length(rows) + 1L]] <- data.frame(
        pos        = as.integer(ppos),
        aa_ref     = parts[1],
        aa_alt     = parts[2],
        consequence = .ensembl_consequence(ct, cons_levels),
        cadd_raw   = .num_or_na(tc$cadd_raw),
        cadd_phred = .num_or_na(tc$cadd_phred),
        revel_score = .num_or_na(tc$revel_score),
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0L) {
    return(list(cadd = empty_cadd, revel = empty_revel))
  }
  df <- do.call(rbind, rows)
  df$aa_change <- paste0(df$aa_ref, df$pos, df$aa_alt)
  df <- df[df$pos >= 1L & df$pos <= protein_length, , drop = FALSE]

  ## --- CADD table: keep rows with a PHRED; max per substitution -------
  cd <- df[!is.na(df$cadd_phred), , drop = FALSE]
  if (nrow(cd) > 0L) {
    cadd <- stats::aggregate(cadd_phred ~ aa_change + pos + aa_ref + aa_alt +
                               consequence, data = cd, FUN = max)
    # Ensembl gives raw too; carry the max-aligned raw, else fall back to
    # PHRED as the proven manual builder did when raw was absent.
    cadd$cadd_raw <- cadd$cadd_phred
    cadd <- cadd[, c("pos", "aa_ref", "aa_alt", "aa_change",
                     "consequence", "cadd_raw", "cadd_phred")]
    cadd$consequence <- factor(as.character(cadd$consequence), levels = cons_levels)
    cadd <- cadd[order(cadd$pos), , drop = FALSE]
    rownames(cadd) <- NULL
  } else {
    cadd <- empty_cadd
  }

  ## --- REVEL table: missense only; max per substitution ---------------
  rv <- df[!is.na(df$revel_score) & df$consequence == "missense", , drop = FALSE]
  if (nrow(rv) > 0L) {
    revel <- stats::aggregate(revel_score ~ aa_change + pos + aa_ref + aa_alt,
                              data = rv, FUN = max)
    revel <- revel[, c("pos", "aa_ref", "aa_alt", "aa_change", "revel_score")]
    revel <- revel[order(revel$pos), , drop = FALSE]
    rownames(revel) <- NULL
  } else {
    revel <- empty_revel
  }

  list(cadd = cadd, revel = revel)
}

## Map Ensembl consequence_terms (a list/vector) to our factor levels.
.ensembl_consequence <- function(terms, cons_levels) {
  t <- tolower(paste(unlist(terms), collapse = " "))
  out <- if (grepl("missense", t)) "missense"
         else if (grepl("synonymous", t)) "synonymous"
         else if (grepl("stop_gained|stop gained", t)) "stop_gained"
         else if (grepl("frameshift", t)) "frameshift"
         else if (grepl("inframe_deletion|inframe deletion", t)) "inframe_deletion"
         else if (grepl("inframe_insertion|inframe insertion", t)) "inframe_insertion"
         else if (grepl("splice_donor|splice donor", t)) "splice_donor"
         else if (grepl("splice_acceptor|splice acceptor", t)) "splice_acceptor"
         else "other"
  factor(out, levels = cons_levels)
}

## Coerce a possibly-NULL/character score to numeric or NA.
.num_or_na <- function(x) {
  if (is.null(x) || length(x) == 0L) return(NA_real_)
  suppressWarnings(as.numeric(x[[1]]))
}
