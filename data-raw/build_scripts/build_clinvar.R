## =================================================================
## build_clinvar.R -- ClinVar variant_summary.txt.gz slicer
## =================================================================

build_clinvar <- function(uniprot_ref, tmp_dir) {

  url <- "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
  local <- file.path(tmp_dir, "variant_summary.txt.gz")
  if (!file.exists(local)) {
    message("  Downloading ClinVar variant_summary.txt.gz (~500 MB) ...")
    utils::download.file(url, local, mode = "wb")
  } else {
    message("  Using cached variant_summary.txt.gz")
  }

  message("  Reading ClinVar ...")
  cv <- readr::read_tsv(local, show_col_types = FALSE,
                         na = c("", "-", "na"),
                         guess_max = 1e5)

  cv <- cv[cv$Assembly == "GRCh38" & !is.na(cv$GeneSymbol), ]
  message(sprintf("  %d GRCh38 entries with gene symbols", nrow(cv)))

  ## ---- Parse HGVS p.* consequence -----------------------------
  three_to_one <- c(Ala="A", Arg="R", Asn="N", Asp="D", Cys="C",
                    Gln="Q", Glu="E", Gly="G", His="H", Ile="I",
                    Leu="L", Lys="K", Met="M", Phe="F", Pro="P",
                    Ser="S", Thr="T", Trp="W", Tyr="Y", Val="V",
                    Ter="*", `=` = "=")

  rx <- "\\(p\\.([A-Za-z]{3})(\\d+)([A-Za-z*]+|fs|del|dup|ins.*?)\\)"
  m <- regmatches(cv$Name, regexec(rx, cv$Name))
  parsed <- do.call(rbind, lapply(m, function(x) {
    if (length(x) < 4L) return(c(NA, NA, NA))
    x[2:4]
  }))

  cv$aa_ref <- unname(three_to_one[parsed[, 1]])
  cv$pos    <- suppressWarnings(as.integer(parsed[, 2]))
  alt_raw   <- parsed[, 3]
  cv$aa_alt <- ifelse(
    alt_raw %in% names(three_to_one),
    unname(three_to_one[alt_raw]),
    tolower(alt_raw)
  )
  cv$aa_change <- ifelse(
    !is.na(cv$aa_ref) & !is.na(cv$pos),
    paste0(cv$aa_ref, cv$pos,
           ifelse(is.na(cv$aa_alt), "?", cv$aa_alt)),
    NA_character_
  )
  # Drop unparseable
  cv <- cv[!is.na(cv$pos) & !is.na(cv$aa_ref), ]

  ## ---- Normalise clinical significance ------------------------
  sig <- cv$ClinicalSignificance
  sig <- gsub("/", "_", sig, fixed = TRUE)
  sig <- gsub(" ", "_", sig, fixed = TRUE)
  sig <- ifelse(grepl("Uncertain_significance", sig), "VUS", sig)
  sig <- ifelse(grepl("Conflicting", sig), "Conflicting", sig)
  sig <- ifelse(sig %in% c("Pathogenic","Likely_pathogenic","VUS",
                            "Likely_benign","Benign","Conflicting"),
                 sig, "Other")
  cv$significance <- factor(sig,
    levels = c("Pathogenic","Likely_pathogenic","VUS",
               "Likely_benign","Benign","Conflicting","Other")
  )

  ## ---- Map ReviewStatus to star factor ------------------------
  star_map <- c(
    "practice guideline"                                    = "4_star",
    "reviewed by expert panel"                              = "3_star",
    "criteria provided, multiple submitters, no conflicts"  = "2_star",
    "criteria provided, single submitter"                   = "1_star",
    "criteria provided, conflicting classifications"        = "1_star",
    "no assertion criteria provided"                        = "0_star",
    "no classification provided"                            = "0_star",
    "no classification for the single variant"              = "0_star"
  )
  cv$review_status <- factor(
    unname(star_map[tolower(cv$ReviewStatus)]),
    levels = c("4_star","3_star","2_star","1_star","0_star")
  )
  cv$review_status[is.na(cv$review_status)] <- "0_star"

  cv$clinvar_id <- paste0("VCV", sprintf("%09d", cv$VariationID))

  ## ---- Per-gene slicing ---------------------------------------
  ## Deduplicate by aa_change + take highest-star entry
  star_pri <- c("4_star"=1, "3_star"=2, "2_star"=3,
                "1_star"=4, "0_star"=5)
  cv$.pri <- star_pri[as.character(cv$review_status)]

  out <- list()
  for (gene in unique(cv$GeneSymbol)) {
    sub <- cv[cv$GeneSymbol == gene,
              c("pos","aa_ref","aa_alt","aa_change",
                "significance","review_status","clinvar_id",".pri")]
    sub <- sub[order(sub$.pri), ]
    sub <- sub[!duplicated(sub$aa_change), ]
    sub$.pri <- NULL
    rownames(sub) <- NULL
    out[[gene]] <- sub
  }
  message(sprintf("  Sliced into %d genes", length(out)))
  out
}
