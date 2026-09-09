## =================================================================
## build_gnomad.R -- gnomAD v4.1 per-variant builder
## =================================================================
##
## Two approaches; choose whichever matches the actual gnomAD v4.1
## release files at build time.
##
## Approach A (preferred): pre-VEP'd coding-variant TSVs from the
## gnomAD downloads page. Stream-parse and slice by gene symbol.
##
## Approach B (fallback): per-gene GraphQL API queries against
## https://gnomad.broadinstitute.org/api. Slow (~6-8 hours for
## 20K genes) but always current. Use this if the bulk files are
## not in a parseable format for the release we want.
##
## This function currently implements Approach B; switch to A by
## populating `tsv_url` and uncommenting the bulk-parse branch.

build_gnomad <- function(uniprot_ref,
                          tmp_dir,
                          sleep_per_query = 0.5,
                          tsv_url = NULL) {

  cons_levels <- c("missense","synonymous","stop_gained",
                   "frameshift","inframe_deletion","inframe_insertion",
                   "splice_donor","splice_acceptor","other")

  ## ---- Approach A: bulk TSV (uncomment when ready) ------------
  ## if (!is.null(tsv_url)) {
  ##   local <- file.path(tmp_dir, basename(tsv_url))
  ##   if (!file.exists(local)) {
  ##     message("  Downloading gnomAD coding TSV ...")
  ##     utils::download.file(tsv_url, local, mode = "wb")
  ##   }
  ##   df <- readr::read_tsv(local, show_col_types = FALSE)
  ##   # ... parse, normalise columns, slice by gene ...
  ## }

  ## ---- Approach B: GraphQL API --------------------------------
  if (!requireNamespace("httr", quietly = TRUE) ||
      !requireNamespace("jsonlite", quietly = TRUE)) {
    stop("build_gnomad needs `httr` and `jsonlite` installed.")
  }
  api <- "https://gnomad.broadinstitute.org/api"
  three_to_one <- c(Ala="A", Arg="R", Asn="N", Asp="D", Cys="C",
                    Gln="Q", Glu="E", Gly="G", His="H", Ile="I",
                    Leu="L", Lys="K", Met="M", Phe="F", Pro="P",
                    Ser="S", Thr="T", Trp="W", Tyr="Y", Val="V",
                    Ter="*")

  query_template <- '{ gene(gene_symbol: "%s", reference_genome: GRCh38) {
      variants(dataset: gnomad_r4) {
        pos consequence hgvsp filters
        exome { ac an af populations { id ac an } }
        genome { ac an af populations { id ac an } }
      }
    } }'

  genes <- unique(uniprot_ref$gene)

  # ---- Per-gene checkpoint cache (makes the crawl resumable) --------
  # Each gene's parsed result is written to <tmp_dir>/gnomad_cache/<gene>.rds
  # as soon as the API returns a definitive 200 (even a 0-row result, so
  # empties are not re-queried). A re-run reads the cache and skips those
  # genes' API calls entirely -- so a crawl that dies at gene 12,000
  # resumes from ~12,000, not from zero. Network failures are NOT cached,
  # so those genes are retried on the next run.
  cache_dir <- file.path(tmp_dir, "gnomad_cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  save_atomic <- function(obj, dest) {
    tmp <- tempfile(tmpdir = dirname(dest), fileext = ".part")
    saveRDS(obj, tmp)
    if (!file.rename(tmp, dest)) {
      file.copy(tmp, dest, overwrite = TRUE)
      unlink(tmp)
    }
  }

  out <- list()
  n_cached <- 0L
  n_fetched <- 0L
  for (i in seq_along(genes)) {
    g <- genes[i]
    if (i %% 200L == 0L)
      message(sprintf("  %d / %d (%s) [fetched %d, from cache %d]",
                      i, length(genes), g, n_fetched, n_cached))

    cache_file <- file.path(cache_dir, paste0(g, ".rds"))
    if (file.exists(cache_file)) {
      # Already fetched on a previous run (df may be 0-row for "no data").
      cached <- tryCatch(readRDS(cache_file), error = function(e) NULL)
      if (is.data.frame(cached) && nrow(cached) > 0L) out[[g]] <- cached
      n_cached <- n_cached + 1L
      next
    }

    # Per-gene tryCatch: a parse/HTTP error on one gene must not halt the
    # whole crawl. The gene's df (parsed below) is assigned via `<<-`.
    df <- NULL
    ok <- tryCatch({
      body <- jsonlite::toJSON(
        list(query = sprintf(query_template, g)),
        auto_unbox = TRUE)
      resp <- httr::POST(api, body = body, encode = "raw",
                         httr::content_type_json(),
                         httr::timeout(60))
      if (httr::status_code(resp) != 200L) {
        # Not a definitive answer (rate-limited/5xx): do NOT cache; retry
        # on the next run.
        FALSE
      } else {
        payload <- httr::content(resp, as = "parsed", simplifyVector = TRUE)
        vs <- payload$data$gene$variants
        df <<- .gnomad_parse_variants(vs, three_to_one, cons_levels)
        TRUE
      }
    }, error = function(e) {
      warning(sprintf("gnomAD query failed for %s: %s",
                      g, conditionMessage(e)))
      FALSE   # transient -> not cached -> retried next run
    })

    if (isTRUE(ok)) {
      # Cache the definitive result (even 0-row) so it is never re-queried.
      save_atomic(df, cache_file)
      n_fetched <- n_fetched + 1L
      if (is.data.frame(df) && nrow(df) > 0L) out[[g]] <- df
    }
    Sys.sleep(sleep_per_query)
  }
  message(sprintf("  Got data for %d genes (%d fetched this run, %d from cache)",
                  length(out), n_fetched, n_cached))
  out
}

# Parse a gnomAD GraphQL `variants` table into the per-gene schema.
# Returns a data.frame (0 rows when `vs` is empty/NULL), never errors on
# an empty result -- so a definitive "no variants" answer is cacheable.
.gnomad_parse_variants <- function(vs, three_to_one, cons_levels) {
    empty <- data.frame(
      pos = integer(), aa_ref = character(), aa_alt = character(),
      aa_change = character(),
      consequence = factor(character(), levels = cons_levels),
      af_exome = numeric(), af_genome = numeric(), af_joint = numeric(),
      ac_joint = integer(), an_joint = integer(), filter = character(),
      stringsAsFactors = FALSE
    )
    if (is.null(vs) || !is.data.frame(vs) || nrow(vs) == 0L) return(empty)

    # Parse HGVSp
    rx <- "^p\\.([A-Z][a-z]{2})(\\d+)([A-Z][a-z]{2}|=|Ter|fs|del|dup|ins.*?)$"
    m <- regmatches(vs$hgvsp, regexec(rx, vs$hgvsp))
    parsed <- do.call(rbind, lapply(m, function(x) {
      if (length(x) < 4L) return(c(NA, NA, NA)); x[2:4]
    }))
    aa_ref <- unname(three_to_one[parsed[, 1]])
    pos    <- suppressWarnings(as.integer(parsed[, 2]))
    aa_alt <- ifelse(parsed[, 3] %in% names(three_to_one),
                     unname(three_to_one[parsed[, 3]]),
                     tolower(parsed[, 3]))
    aa_change <- ifelse(!is.na(aa_ref) & !is.na(pos),
                        paste0(aa_ref, pos,
                               ifelse(is.na(aa_alt), "?", aa_alt)),
                        NA_character_)

    cons <- vs$consequence
    cons_norm <- ifelse(cons %in% cons_levels, cons, "other")

    # Joint AF -- use exome if available, else genome
    af_ex <- if (is.data.frame(vs$exome)) vs$exome$af else NA_real_
    af_gn <- if (is.data.frame(vs$genome)) vs$genome$af else NA_real_
    ac_ex <- if (is.data.frame(vs$exome)) vs$exome$ac else NA_integer_
    an_ex <- if (is.data.frame(vs$exome)) vs$exome$an else NA_integer_
    ac_gn <- if (is.data.frame(vs$genome)) vs$genome$ac else NA_integer_
    an_gn <- if (is.data.frame(vs$genome)) vs$genome$an else NA_integer_

    # joint = exome + genome counts where both present
    ac_joint <- ifelse(is.na(ac_ex), 0L, ac_ex) +
                 ifelse(is.na(ac_gn), 0L, ac_gn)
    an_joint <- ifelse(is.na(an_ex), 0L, an_ex) +
                 ifelse(is.na(an_gn), 0L, an_gn)
    af_joint <- ifelse(an_joint > 0L, ac_joint / an_joint, NA_real_)
    ac_joint <- as.integer(ac_joint)
    an_joint <- as.integer(an_joint)

    filt <- vs$filters
    filt <- ifelse(vapply(filt, length, integer(1)) == 0L, "PASS",
                    vapply(filt, paste, character(1), collapse = ","))

    df <- data.frame(
      pos = pos,
      aa_ref = aa_ref,
      aa_alt = aa_alt,
      aa_change = aa_change,
      consequence = factor(cons_norm, levels = cons_levels),
      af_exome = af_ex,
      af_genome = af_gn,
      af_joint = af_joint,
      ac_joint = ac_joint,
      an_joint = an_joint,
      filter = filt,
      stringsAsFactors = FALSE
    )
    df <- df[!is.na(df$pos) & !is.na(df$aa_change), ]
    rownames(df) <- NULL
    df
}
