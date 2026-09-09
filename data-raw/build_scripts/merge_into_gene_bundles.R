## =================================================================
## merge_into_gene_bundles.R
## =================================================================
##
## Takes the six per-source intermediates (each a named list keyed
## by HGNC symbol) and produces one combined per-gene .rds file
## conforming to DATA_FORMAT_SPEC.md. Validates each bundle before
## writing.

merge_into_gene_bundles <- function(out_dir,
                                     uniprot_ref,
                                     domains,
                                     clinvar,
                                     gnomad,
                                     alphamissense,
                                     revel,
                                     cadd) {
  build_date <- Sys.Date()
  cons_levels <- c("missense","synonymous","stop_gained","frameshift",
                   "inframe_deletion","inframe_insertion",
                   "splice_donor","splice_acceptor","other")
  sig_levels  <- c("Pathogenic","Likely_pathogenic","VUS",
                   "Likely_benign","Benign","Conflicting","Other")
  star_levels <- c("4_star","3_star","2_star","1_star","0_star")
  src_levels  <- c("Pfam","InterPro","SMART","PROSITE","PRINTS","PANTHER")
  am_levels   <- c("likely_benign","ambiguous","likely_pathogenic")

  # Helper: an empty data.frame conforming to a slot's schema
  empty_df <- function(slot) {
    switch(slot,
      "domains" = data.frame(
        start = integer(), end = integer(),
        name = character(), accession = character(),
        source = factor(character(), levels = src_levels),
        stringsAsFactors = FALSE
      ),
      "clinvar" = data.frame(
        pos = integer(), aa_ref = character(), aa_alt = character(),
        aa_change = character(),
        significance = factor(character(), levels = sig_levels),
        review_status = factor(character(), levels = star_levels),
        clinvar_id = character(),
        stringsAsFactors = FALSE
      ),
      "gnomad" = data.frame(
        pos = integer(), aa_ref = character(), aa_alt = character(),
        aa_change = character(),
        consequence = factor(character(), levels = cons_levels),
        af_joint = numeric(), ac_joint = integer(),
        an_joint = integer(), filter = character(),
        stringsAsFactors = FALSE
      ),
      "alphamissense" = data.frame(
        pos = integer(), aa_ref = character(), aa_alt = character(),
        aa_change = character(), am_score = numeric(),
        am_class = factor(character(), levels = am_levels),
        stringsAsFactors = FALSE
      ),
      "revel" = data.frame(
        pos = integer(), aa_ref = character(), aa_alt = character(),
        aa_change = character(), revel_score = numeric(),
        stringsAsFactors = FALSE
      ),
      "cadd" = data.frame(
        pos = integer(), aa_ref = character(), aa_alt = character(),
        aa_change = character(),
        consequence = factor(character(), levels = cons_levels),
        cadd_raw = numeric(), cadd_phred = numeric(),
        stringsAsFactors = FALSE
      )
    )
  }

  # Source-version string (filled by individual builders; if not
  # present, default to date)
  source_versions <- sprintf(
    "InterPro=%s;ClinVar=%s;gnomAD=v4.1;AlphaMissense=2023-09;REVEL=v1.3;CADD=v1.7",
    build_date, build_date)

  n_written <- 0L
  n_skipped <- 0L    # already present from a previous run (resume)
  n_failed  <- 0L    # validation failure or a hard error in the gene body
  failed_genes <- character()

  # Atomic write: saveRDS to a temp file in the same dir, then rename, so
  # a crash mid-write never leaves a truncated <gene>.rds that the resume
  # skip-guard would mistake for a completed gene.
  save_atomic <- function(obj, dest) {
    tmp <- tempfile(tmpdir = dirname(dest), fileext = ".rds.part")
    saveRDS(obj, tmp, compress = "xz")
    if (!file.rename(tmp, dest)) {
      file.copy(tmp, dest, overwrite = TRUE)
      unlink(tmp)
    }
  }

  # Build (and write) one gene's bundle. Returns "written", "skipped",
  # or "failed". Any error inside is caught by the caller's tryCatch.
  build_one <- function(i) {
    row <- uniprot_ref[i, ]
    gene <- row$gene
    dest <- file.path(out_dir, paste0(gene, ".rds"))

    # Resume: a completed gene from a prior run is left untouched.
    if (file.exists(dest)) {
      return("skipped")
    }

    meta <- data.frame(
      gene                  = gene,
      uniprot_id            = row$uniprot_id,
      protein_length        = as.integer(row$protein_length),
      ensembl_gene_id       = if (is.na(row$ensembl_gene_id)) ""
                              else row$ensembl_gene_id,
      ensembl_transcript_id = if (is.na(row$ensembl_transcript_id)) ""
                              else row$ensembl_transcript_id,
      build_date            = build_date,
      source_versions       = source_versions,
      stringsAsFactors      = FALSE
    )

    bundle <- list(
      meta          = meta,
      domains       = if (!is.null(domains[[gene]])) domains[[gene]]
                      else empty_df("domains"),
      clinvar       = if (!is.null(clinvar[[gene]])) clinvar[[gene]]
                      else empty_df("clinvar"),
      gnomad        = if (!is.null(gnomad[[gene]])) gnomad[[gene]]
                      else empty_df("gnomad"),
      alphamissense = if (!is.null(alphamissense[[gene]])) alphamissense[[gene]]
                      else empty_df("alphamissense"),
      revel         = if (!is.null(revel[[gene]])) revel[[gene]]
                      else empty_df("revel"),
      cadd          = if (!is.null(cadd[[gene]])) cadd[[gene]]
                      else empty_df("cadd")
    )

    v <- validate_gene_data(bundle, strict = FALSE)
    if (!v$valid) {
      warning(sprintf("Validation failed for %s; skipping.\n  %s",
                      gene, paste(v$issues, collapse = "\n  ")))
      return("failed")
    }

    save_atomic(bundle, dest)
    "written"
  }

  for (i in seq_len(nrow(uniprot_ref))) {
    gene <- uniprot_ref$gene[i]
    if (i %% 500L == 0L)
      message(sprintf("    [%d / %d] %s", i, nrow(uniprot_ref), gene))

    # Per-gene tryCatch: a hard error (bad parse, I/O, etc.) on ONE gene
    # is logged and skipped, never halting a multi-hour whole-genome run.
    status <- tryCatch(
      build_one(i),
      error = function(e) {
        warning(sprintf("Error building %s; skipping. %s",
                        gene, conditionMessage(e)))
        "failed"
      }
    )

    if (status == "written") {
      n_written <- n_written + 1L
    } else if (status == "skipped") {
      n_skipped <- n_skipped + 1L
    } else {
      n_failed <- n_failed + 1L
      failed_genes <- c(failed_genes, gene)
    }
  }

  message(sprintf(
    "\n  Wrote %d gene files; %d skipped (already present); %d failed",
    n_written, n_skipped, n_failed))
  if (n_failed > 0L) {
    writeLines(failed_genes, file.path(out_dir, "failed_genes.txt"))
    message("  Failed gene list: ", file.path(out_dir, "failed_genes.txt"))
  }
  invisible(list(written = n_written, skipped = n_skipped,
                 failed = n_failed, failed_genes = failed_genes))
}
