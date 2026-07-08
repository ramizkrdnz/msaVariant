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
  out <- list()
  for (i in seq_along(genes)) {
    g <- genes[i]
    if (i %% 200L == 0L)
      message(sprintf("  %d / %d (%s)", i, length(genes), g))

    body <- jsonlite::toJSON(
      list(query = sprintf(query_template, g)),
      auto_unbox = TRUE)
    resp <- tryCatch(
      httr::POST(api, body = body, encode = "raw",
                 httr::content_type_json(),
                 httr::timeout(60)),
      error = function(e) NULL)
    if (is.null(resp) || httr::status_code(resp) != 200L) {
      Sys.sleep(sleep_per_query); next
    }
    payload <- httr::content(resp, as = "parsed",
                              simplifyVector = TRUE)
    vs <- payload$data$gene$variants
    if (is.null(vs) || nrow(vs) == 0L) {
      Sys.sleep(sleep_per_query); next
    }

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
    if (nrow(df) > 0L) out[[g]] <- df
    Sys.sleep(sleep_per_query)
  }
  message(sprintf("  Got data for %d genes", length(out)))
  out
}
