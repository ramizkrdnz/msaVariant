## =================================================================
## msaVariant Zenodo deposit builder (combined per-gene files)
## =================================================================
##
## Produces ~19,000 .rds files (one per HGNC-approved human gene),
## each containing the full per-gene annotation bundle per
## DATA_FORMAT_SPEC.md.
##
## RUN ON: a networked machine with ~500 GB free disk, R ≥ 4.1, and
##   ~6-12 hours of mostly-unattended time.
##
## OUTPUT: zenodo_payload/ containing one <GENE>.rds per gene plus
##   MANIFEST.tsv, README.md, LICENSES/.
##
## STRATEGY:
##   1. Build a UniProt-keyed reference of all human SwissProt
##      proteins (gene symbol, length, transcript IDs).
##   2. For each annotation source, build a per-source intermediate
##      table keyed by HGNC gene symbol.
##   3. Merge the six intermediates into per-gene bundles, validate
##      each against the spec, and write to .rds.
##   4. Package licenses and manifest.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

# Source the per-source builders
SCRIPT_DIR <- file.path("data-raw", "build_scripts")
for (f in c("build_uniprot_ref.R",
            "build_domains.R",
            "build_clinvar.R",
            "build_gnomad.R",
            "build_alphamissense.R",
            "build_revel.R",
            "build_cadd.R",
            "merge_into_gene_bundles.R",
            "build_manifest.R")) {
  source(file.path(SCRIPT_DIR, f))
}

# Also need the validator from the package itself
source(file.path("R", "validate_gene_data.R"))

OUT_DIR <- "zenodo_payload"
TMP_DIR <- "build_tmp"
INTERMEDIATE_DIR <- file.path(TMP_DIR, "intermediates")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TMP_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(INTERMEDIATE_DIR, showWarnings = FALSE, recursive = TRUE)

t0 <- Sys.time()

message("\n=== msaVariant Zenodo deposit build ===\n")
message("Start time: ", format(t0), "\n")

# 1. UniProt reference -- everything else hangs on this
message("Step 1/9: UniProt human SwissProt reference ...")
uniprot_ref <- build_uniprot_ref(tmp_dir = TMP_DIR)
saveRDS(uniprot_ref, file.path(INTERMEDIATE_DIR, "uniprot_ref.rds"))
message(sprintf("  -> %d human SwissProt entries\n", nrow(uniprot_ref)))

# 2-7. Build per-source intermediates
message("Step 2/9: InterPro/Pfam domains ...")
domains_by_gene <- build_domains(uniprot_ref = uniprot_ref,
                                  tmp_dir = TMP_DIR)
saveRDS(domains_by_gene, file.path(INTERMEDIATE_DIR, "domains.rds"))

message("\nStep 3/9: ClinVar ...")
clinvar_by_gene <- build_clinvar(uniprot_ref = uniprot_ref,
                                  tmp_dir = TMP_DIR)
saveRDS(clinvar_by_gene, file.path(INTERMEDIATE_DIR, "clinvar.rds"))

message("\nStep 4/9: gnomAD v4.1 ...")
gnomad_by_gene <- build_gnomad(uniprot_ref = uniprot_ref,
                                tmp_dir = TMP_DIR)
saveRDS(gnomad_by_gene, file.path(INTERMEDIATE_DIR, "gnomad.rds"))

message("\nStep 5/9: AlphaMissense ...")
alphamissense_by_gene <- build_alphamissense(uniprot_ref = uniprot_ref,
                                              tmp_dir = TMP_DIR)
saveRDS(alphamissense_by_gene, file.path(INTERMEDIATE_DIR, "alphamissense.rds"))

message("\nStep 6/9: REVEL ...")
revel_by_gene <- build_revel(uniprot_ref = uniprot_ref,
                              tmp_dir = TMP_DIR)
saveRDS(revel_by_gene, file.path(INTERMEDIATE_DIR, "revel.rds"))

message("\nStep 7/9: CADD ...")
cadd_by_gene <- build_cadd(uniprot_ref = uniprot_ref,
                            tmp_dir = TMP_DIR)
saveRDS(cadd_by_gene, file.path(INTERMEDIATE_DIR, "cadd.rds"))

# 8. Merge into per-gene bundles and validate
message("\nStep 8/9: merging into per-gene bundles ...")
merge_into_gene_bundles(
  out_dir       = OUT_DIR,
  uniprot_ref   = uniprot_ref,
  domains       = domains_by_gene,
  clinvar       = clinvar_by_gene,
  gnomad        = gnomad_by_gene,
  alphamissense = alphamissense_by_gene,
  revel         = revel_by_gene,
  cadd          = cadd_by_gene
)

# 9. Manifest + licenses
message("\nStep 9/9: manifest, licenses, README ...")
build_manifest(out_dir = OUT_DIR)

t1 <- Sys.time()
message("\n=== Build complete ===")
message("End time: ", format(t1))
message("Total elapsed: ",
        format(round(difftime(t1, t0, units = "mins"), 1)))
message("\nNext step: follow data-raw/ZENODO_UPLOAD.md to upload.")
