## =================================================================
## msaVariant: Zenodo deposit builder (v0.1, combined per-gene)
## =================================================================
##
## RUN ON A NETWORKED MACHINE.
##
## Produces one .rds per HGNC gene, conforming to DATA_FORMAT_SPEC.md.
## Each file contains all 6 annotation types plus metadata.
##
## Time required:  6-12 hours
## Disk required:  ~500 GB during build, ~10 GB final deposit
## RAM required:   ~16 GB
##
## Output:
##   zenodo_payload/
##     PATL1.rds
##     TP53.rds
##     ...
##     MANIFEST.tsv
##     README.md
##     LICENSES/

suppressPackageStartupMessages({
  library(readr)
})

OUT_DIR <- "zenodo_payload"
TMP_DIR <- "build_tmp"
SCRIPT_DIR <- file.path("data-raw", "build_scripts")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TMP_DIR, showWarnings = FALSE, recursive = TRUE)

cat("msaVariant Zenodo deposit builder (v0.1)\n")
cat("==========================================\n\n")

## Source the package's own validator -- the merger needs it
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)

## Step 1: UniProt human reference (gene <-> uniprot <-> ensembl)
cat("Step 1/8: UniProt human reference ...\n")
source(file.path(SCRIPT_DIR, "build_uniprot_ref.R"))
uniprot_ref <- build_uniprot_ref(tmp_dir = TMP_DIR)
saveRDS(uniprot_ref, file.path(TMP_DIR, "uniprot_ref.rds"))

## Step 2: InterPro/Pfam domains, sliced by gene
cat("\nStep 2/8: InterPro/Pfam domains ...\n")
source(file.path(SCRIPT_DIR, "build_domains.R"))
domains_by_gene <- build_domains(tmp_dir = TMP_DIR,
                                  uniprot_ref = uniprot_ref)
saveRDS(domains_by_gene, file.path(TMP_DIR, "domains_by_gene.rds"))

## Step 3: ClinVar
cat("\nStep 3/8: ClinVar ...\n")
source(file.path(SCRIPT_DIR, "build_clinvar.R"))
clinvar_by_gene <- build_clinvar(tmp_dir = TMP_DIR,
                                  uniprot_ref = uniprot_ref)
saveRDS(clinvar_by_gene, file.path(TMP_DIR, "clinvar_by_gene.rds"))

## Step 4: gnomAD v4.1
cat("\nStep 4/8: gnomAD v4.1 ...\n")
source(file.path(SCRIPT_DIR, "build_gnomad.R"))
gnomad_by_gene <- build_gnomad(tmp_dir = TMP_DIR,
                                uniprot_ref = uniprot_ref)
saveRDS(gnomad_by_gene, file.path(TMP_DIR, "gnomad_by_gene.rds"))

## Step 5: AlphaMissense
cat("\nStep 5/8: AlphaMissense ...\n")
source(file.path(SCRIPT_DIR, "build_alphamissense.R"))
am_by_gene <- build_alphamissense(tmp_dir = TMP_DIR,
                                   uniprot_ref = uniprot_ref)
saveRDS(am_by_gene, file.path(TMP_DIR, "alphamissense_by_gene.rds"))

## Step 6: REVEL
cat("\nStep 6/8: REVEL ...\n")
source(file.path(SCRIPT_DIR, "build_revel.R"))
revel_by_gene <- build_revel(tmp_dir = TMP_DIR,
                              uniprot_ref = uniprot_ref)
saveRDS(revel_by_gene, file.path(TMP_DIR, "revel_by_gene.rds"))

## Step 7: CADD v1.7
cat("\nStep 7/8: CADD v1.7 ...\n")
source(file.path(SCRIPT_DIR, "build_cadd.R"))
cadd_by_gene <- build_cadd(tmp_dir = TMP_DIR,
                            uniprot_ref = uniprot_ref)
saveRDS(cadd_by_gene, file.path(TMP_DIR, "cadd_by_gene.rds"))

## Step 8: merge into combined per-gene bundles, validate, write
cat("\nStep 8/8: merging into per-gene bundles ...\n")
source(file.path(SCRIPT_DIR, "merge_into_gene_bundles.R"))
merge_into_gene_bundles(
  out_dir       = OUT_DIR,
  uniprot_ref   = uniprot_ref,
  domains       = domains_by_gene,
  clinvar       = clinvar_by_gene,
  gnomad        = gnomad_by_gene,
  alphamissense = am_by_gene,
  revel         = revel_by_gene,
  cadd          = cadd_by_gene
)

## Generate the manifest
cat("\nGenerating MANIFEST ...\n")
source(file.path(SCRIPT_DIR, "build_manifest.R"))
build_manifest(payload_dir = OUT_DIR)

cat("\nDone. Inspect '", OUT_DIR, "/', then follow ZENODO_UPLOAD.md.\n",
    sep = "")
