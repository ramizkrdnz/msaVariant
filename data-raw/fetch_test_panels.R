## =================================================================
## msaVariant test cases — three independent real MSAs
## =================================================================
##
## This script does NOT do the alignment for you. It downloads the
## unaligned ortholog sequences for three test genes from UniProt
## and saves them to FASTA files. You then run alignment + plotting
## yourself, following the workflow in `tutorial_full_workflow.Rmd`.
##
## Why we don't ship pre-aligned sequences in this script:
##   The alignment is part of what you're testing. If we did it
##   ahead of time you'd be testing the visualisation, not the full
##   pipeline.
##
## Three independent genes are included so any quirks of a single
## protein don't mask bugs:
##
##   1. TP53     — small (393 aa), heavily annotated tumor-suppressor.
##                 Variant of interest: p.R175H (a recurrent hotspot).
##   2. CFTR     — large (1480 aa), common CF variant p.F508del.
##                 Tests indel handling and large-alignment windows.
##   3. BRCA1    — very large (1863 aa) with BRCT domains and many
##                 VUS in ClinVar. Tests the canonical use case.
##
## Each gene is fetched across the same eight species so you can
## compare the visualisations directly.

suppressPackageStartupMessages({
  library(Biostrings)
})

## --- Helper: pull one UniProt sequence by accession ---------------
## Returns an AAStringSet of length 1. Errors clearly if the
## accession doesn't resolve or the canonical sequence is empty.
fetch_uniprot <- function(accession) {
  url <- sprintf("https://rest.uniprot.org/uniprotkb/%s.fasta", accession)
  tryCatch({
    aa <- readAAStringSet(url)
    if (length(aa) == 0L) {
      stop("UniProt returned no sequence for ", accession)
    }
    aa[1L]   # canonical
  }, error = function(e) {
    message("FAILED to fetch ", accession, ": ", conditionMessage(e))
    NULL
  })
}

## --- Helper: fetch a panel of orthologs --------------------------
## `accs` is a named character vector. Names become FASTA headers.
fetch_panel <- function(accs, out_path) {
  seqs <- lapply(accs, fetch_uniprot)

  # Validate: each element must be an AAStringSet of length >= 1.
  # Anything else (NULL from a failed fetch, or an unexpected
  # return type) gets dropped with a clear message.
  is_valid <- vapply(seqs, function(x) {
    inherits(x, "AAStringSet") && length(x) >= 1L
  }, logical(1L))

  if (any(!is_valid)) {
    message("Dropped ", sum(!is_valid),
            " invalid/empty sequence(s): ",
            paste(names(accs)[!is_valid], collapse = ", "))
  }
  if (sum(is_valid) < 2L) {
    message("Not enough valid sequences fetched; skipping ", out_path)
    return(invisible(NULL))
  }

  # Use Reduce(append, ...) for safe AAStringSet concatenation.
  # `do.call(c, ...)` can silently degrade to a plain list when a
  # list element isn't the expected class, which is exactly the
  # bug this avoids.
  valid_seqs <- seqs[is_valid]
  out <- Reduce(append, valid_seqs)

  # Apply readable names from the input
  names(out) <- names(accs)[is_valid]

  writeXStringSet(out, out_path)
  message(sprintf("Wrote %d sequences to %s", length(out), out_path))
  invisible(out)
}

## =================================================================
## TEST CASE 1: TP53
## =================================================================
## Tumor protein p53. 393 aa in human. Most-studied tumor-suppressor
## with thousands of ClinVar entries. Hotspot variants include
## R175H, R248Q/W, R273H, R282W — all in the DNA-binding domain
## (residues ~94-312).
tp53 <- c(
  TP53_HUMAN = "P04637",
  TP53_PANTR = "P61260",
  TP53_MOUSE = "P02340",
  TP53_RAT   = "P10361",
  TP53_BOVIN = "Q9TUB2",
  TP53_CHICK = "P10360",
  TP53_XENLA = "P07193",
  TP53_DANRE = "P79734"
)
fetch_panel(tp53, "tp53_unaligned.fasta")

## =================================================================
## TEST CASE 2: CFTR
## =================================================================
## Cystic fibrosis transmembrane conductance regulator. 1480 aa in
## human. p.F508del (a 3-nt deletion removing F508) is the most
## common CF-causing variant globally. Tests indel handling in
## the variant overlay and long-alignment windowing.
cftr <- c(
  CFTR_HUMAN = "P13569",
  CFTR_PANTR = "Q4H4M2",
  CFTR_MOUSE = "P26361",
  CFTR_RAT   = "P34158",
  CFTR_BOVIN = "P35071",
  CFTR_CHICK = "Q98SX8",
  CFTR_XENLA = "P26361",         # NOTE: Xenopus CFTR is often
                                  # incomplete in UniProt; the fetch
                                  # may silently skip this entry.
                                  # That's fine — the panel still
                                  # has 7+ species.
  CFTR_DANRE = "Q9I8G0"
)
fetch_panel(cftr, "cftr_unaligned.fasta")

## =================================================================
## TEST CASE 3: BRCA1
## =================================================================
## Breast cancer type 1 susceptibility protein. 1863 aa in human.
## RING domain (~1-100) and BRCT repeats (~1646-1859). High VUS
## burden makes it a real-world clinical-interpretation testbed.
brca1 <- c(
  BRCA1_HUMAN = "P38398",
  BRCA1_PANTR = "Q9GKK8",
  BRCA1_MOUSE = "P48754",
  BRCA1_RAT   = "O54952",
  BRCA1_BOVIN = "Q865B0",
  BRCA1_CHICK = "Q90842",
  BRCA1_XENLA = "Q6PFM4",
  BRCA1_DANRE = "Q5RHM4"
)
fetch_panel(brca1, "brca1_unaligned.fasta")

message("\nDONE. Three unaligned FASTA files written:")
message("  tp53_unaligned.fasta")
message("  cftr_unaligned.fasta")
message("  brca1_unaligned.fasta")
message("\nNow run alignment + plotting as shown in ",
        "vignettes/tutorial_full_workflow.Rmd")
message("\nSuggested variants of interest for each:")
message("  TP53  p.R175H   (hotspot in DNA-binding domain)")
message("  CFTR  p.F508del (most common CF variant)")
message("  BRCA1 p.C61G    (RING-domain pathogenic VUS-resolved)")
