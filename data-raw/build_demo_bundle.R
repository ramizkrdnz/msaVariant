## =====================================================================
## build_demo_bundle.R
## ---------------------------------------------------------------------
## Builds a SYNTHETIC, license-clean example gene bundle ("DEMO1") plus a
## matching synthetic ortholog alignment, and ships them in inst/extdata/
## so the package's Quick Example / man examples run out-of-the-box with
## NO Zenodo download and NO real (license-encumbered) annotation data.
##
## Everything here is fabricated. The numbers are illustrative only and
## must never be read as real AlphaMissense / CADD / REVEL / ClinVar /
## gnomAD predictions. DEMO1 is not an HGNC gene symbol.
##
## The synthetic variant DEMO1 p.R21H is engineered so that all five
## ACMG codes fire (PS1, PM1, PM2, PM5, PP3), reproducing the full strip.
##
## Run from the package root with:  Rscript data-raw/build_demo_bundle.R
## =====================================================================

suppressMessages(devtools::load_all(".", quiet = TRUE))
set.seed(42)

## --- Canonical factor levels (must match validate_gene_data() exactly)
SIGNIFICANCE_LEVELS <- c("Pathogenic", "Likely_pathogenic", "VUS",
                         "Likely_benign", "Benign", "Conflicting", "Other")
REVIEW_STATUS_LEVELS <- c("4_star", "3_star", "2_star", "1_star", "0_star")
DOMAIN_SOURCE_LEVELS <- c("Pfam", "InterPro", "SMART", "PROSITE",
                          "PRINTS", "PANTHER")
CONSEQUENCE_LEVELS <- c("missense", "synonymous", "stop_gained",
                        "frameshift", "inframe_deletion",
                        "inframe_insertion", "splice_donor",
                        "splice_acceptor", "other")
AM_CLASS_LEVELS <- c("likely_benign", "ambiguous", "likely_pathogenic")

AA <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]

## --- Reference protein (40 aa, gapless; position 21 == 'R') ----------
REF <- "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRA"
L <- nchar(REF)                       # 40
stopifnot(L == 40L, substr(REF, 21, 21) == "R")
ref_aa <- function(pos) substr(REF, pos, pos)
alt_for <- function(pos) {            # a deterministic "alt" != ref
  r <- ref_aa(pos); a <- if (r != "A") "A" else "G"; a
}
aac <- function(pos, alt = NULL) {
  alt <- if (is.null(alt)) alt_for(pos) else alt
  sprintf("%s%d%s", ref_aa(pos), pos, alt)
}

VPOS <- 21L                            # the showcase variant position
VALT <- "H"                            # DEMO1 p.R21H

## --- meta (1 row) ----------------------------------------------------
meta <- data.frame(
  gene                  = "DEMO1",
  uniprot_id            = "DEMO00001",
  protein_length        = L,
  ensembl_gene_id       = "ENSGDEMO00000000001",
  ensembl_transcript_id = "ENSTDEMO00000000001",
  refseq_protein        = "NP_DEMO0001.1",
  mane_select           = "NM_DEMO0001.1",
  build_date            = Sys.Date(),
  source_versions       = paste0(
    "SYNTHETIC — illustrative only, not real ",
    "AlphaMissense/CADD/ClinVar/gnomAD"),
  stringsAsFactors = FALSE
)

## --- domains (variant sits inside one -> PM1) ------------------------
domains <- data.frame(
  start     = c(5L, 2L),
  end       = c(35L, 10L),
  name      = c("Demo functional domain", "Demo N-terminal region"),
  accession = c("PFDEMO0001", "IPRDEMO001"),
  source    = factor(c("Pfam", "InterPro"), levels = DOMAIN_SOURCE_LEVELS),
  stringsAsFactors = FALSE
)

## --- clinvar (~30 rows; includes PS1 + PM5 drivers) ------------------
cv_pos <- sort(unique(c(VPOS, sample(setdiff(1:L, VPOS), 26))))
clinvar <- do.call(rbind, lapply(seq_along(cv_pos), function(i) {
  p <- cv_pos[i]
  sig <- sample(SIGNIFICANCE_LEVELS, 1)
  data.frame(pos = as.integer(p), aa_ref = ref_aa(p), aa_alt = alt_for(p),
             aa_change = aac(p), significance = sig,
             review_status = sample(REVIEW_STATUS_LEVELS, 1),
             clinvar_id = sprintf("DEMO%07d", i),
             stringsAsFactors = FALSE)
}))
## PS1 driver: exact aa_change R21H, Pathogenic
clinvar <- rbind(clinvar, data.frame(
  pos = VPOS, aa_ref = "R", aa_alt = VALT, aa_change = sprintf("R%dH", VPOS),
  significance = "Pathogenic", review_status = "2_star",
  clinvar_id = "DEMO9000001", stringsAsFactors = FALSE))
## PM5 driver: same pos, DIFFERENT alt, (Likely_)pathogenic
clinvar <- rbind(clinvar, data.frame(
  pos = VPOS, aa_ref = "R", aa_alt = "C", aa_change = sprintf("R%dC", VPOS),
  significance = "Likely_pathogenic", review_status = "1_star",
  clinvar_id = "DEMO9000002", stringsAsFactors = FALSE))
clinvar$significance  <- factor(clinvar$significance,  levels = SIGNIFICANCE_LEVELS)
clinvar$review_status <- factor(clinvar$review_status, levels = REVIEW_STATUS_LEVELS)
rownames(clinvar) <- NULL

## --- gnomad (~24 rows; DELIBERATELY no R21H -> PM2 as "absent") ------
gn_pos <- sort(sample(setdiff(1:L, VPOS), 24))
gnomad <- data.frame(
  pos = as.integer(gn_pos),
  aa_ref = vapply(gn_pos, ref_aa, ""),
  aa_alt = vapply(gn_pos, alt_for, ""),
  aa_change = vapply(gn_pos, aac, ""),
  consequence = factor(sample(c("missense", "synonymous"), length(gn_pos),
                              replace = TRUE, prob = c(0.8, 0.2)),
                       levels = CONSEQUENCE_LEVELS),
  af_joint = round(runif(length(gn_pos), 1e-6, 5e-3), 7),
  ac_joint = as.integer(sample(1:40, length(gn_pos), replace = TRUE)),
  an_joint = rep(152312L, length(gn_pos)),
  filter = "PASS",
  stringsAsFactors = FALSE
)

## --- alphamissense / revel / cadd: one row per position -------------
am_score <- runif(L, 0, 1); am_score[VPOS] <- 0.98
am_class <- ifelse(am_score > 0.564, "likely_pathogenic",
                   ifelse(am_score < 0.34, "likely_benign", "ambiguous"))
am_class[VPOS] <- "likely_pathogenic"
alphamissense <- data.frame(
  pos = 1:L, aa_ref = vapply(1:L, ref_aa, ""),
  aa_alt = c(vapply(1:(VPOS-1), alt_for, ""), VALT, vapply((VPOS+1):L, alt_for, "")),
  aa_change = c(vapply(1:(VPOS-1), aac, ""), sprintf("R%dH", VPOS),
                vapply((VPOS+1):L, aac, "")),
  am_score = round(am_score, 4),
  am_class = factor(am_class, levels = AM_CLASS_LEVELS),
  stringsAsFactors = FALSE
)
revel_score <- runif(L, 0, 1); revel_score[VPOS] <- 0.92
revel <- data.frame(
  pos = 1:L, aa_ref = alphamissense$aa_ref, aa_alt = alphamissense$aa_alt,
  aa_change = alphamissense$aa_change, revel_score = round(revel_score, 4),
  stringsAsFactors = FALSE)
cadd_phred <- runif(L, 0, 40); cadd_phred[VPOS] <- 32
cadd <- data.frame(
  pos = 1:L, aa_ref = alphamissense$aa_ref, aa_alt = alphamissense$aa_alt,
  aa_change = alphamissense$aa_change,
  consequence = factor(rep("missense", L), levels = CONSEQUENCE_LEVELS),
  cadd_raw = round(cadd_phred / 8, 3), cadd_phred = round(cadd_phred, 2),
  stringsAsFactors = FALSE)
alphamissense$pos <- as.integer(alphamissense$pos)
revel$pos <- as.integer(revel$pos)
cadd$pos  <- as.integer(cadd$pos)

## --- Assemble + validate --------------------------------------------
bundle <- list(meta = meta, domains = domains, clinvar = clinvar,
               gnomad = gnomad, alphamissense = alphamissense,
               revel = revel, cadd = cadd)
v <- validate_gene_data(bundle, strict = FALSE)
if (!isTRUE(v$valid)) stop("DEMO1 bundle FAILED validation:\n  ",
                           paste(v$issues, collapse = "\n  "))
message("DEMO1 bundle validates OK.")
codes <- compute_acmg_codes(bundle, VPOS, sprintf("R%dH", VPOS))
message("ACMG codes for DEMO1 p.R", VPOS, "H: ", paste(codes, collapse = ", "))
stopifnot(identical(codes, c("PS1", "PM1", "PM2", "PM5", "PP3")))

dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
saveRDS(bundle, "inst/extdata/DEMO1.rds")
message("Wrote inst/extdata/DEMO1.rds")

## --- Synthetic ortholog alignment (8 species, 40 cols, gapless ref) --
species <- c("DEMO1_HUMAN", "DEMO1_PANTR", "DEMO1_MOUSE", "DEMO1_RAT",
             "DEMO1_BOVIN", "DEMO1_CHICK", "DEMO1_XENLA", "DEMO1_DANRE")
seqs <- setNames(vector("list", length(species)), species)
seqs[["DEMO1_HUMAN"]] <- strsplit(REF, "")[[1]]
for (sp in species[-1]) {
  s <- strsplit(REF, "")[[1]]
  muts <- sample(setdiff(1:L, VPOS), sample(5:9, 1))   # keep pos 21 conserved
  for (m in muts) s[m] <- sample(setdiff(AA, s[m]), 1)
  if (runif(1) < 0.4) s[sample(setdiff(1:L, VPOS), 1)] <- "-"  # occasional gap
  seqs[[sp]] <- s
}
fa <- unlist(lapply(species, function(sp)
  c(paste0(">", sp), paste(seqs[[sp]], collapse = ""))))
writeLines(fa, "inst/extdata/demo_aligned.fasta")
message("Wrote inst/extdata/demo_aligned.fasta")
message("Done.")
