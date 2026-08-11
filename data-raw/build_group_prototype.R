## =================================================================
## build_group_prototype.R
## =================================================================
##
## Prototype of the Path 2 "grouped-bundle" architecture at small
## scale. Packs multiple per-gene bundles into a few group .rds files,
## writes the gene->group index that ships in the package, and builds
## a MANIFEST over the group files.
##
## This is a TEST-SCALE prototype (5 genes, 2 groups), not the real
## ~19,500-gene build. Run from the package root:
##   Rscript data-raw/build_group_prototype.R
##
## Inputs (already in the repo):
##   tests/testthat/fixtures/TP53.rds, WDR31.rds   (real bundles)
##   inst/extdata/DEMO1.rds                        (synthetic bundle)
## Outputs:
##   data-raw/zenodo_grouped_payload/group_01.rds  {TP53, WDR31}
##   data-raw/zenodo_grouped_payload/group_02.rds  {DEMO1, SYNTH1, SYNTH2}
##   data-raw/zenodo_grouped_payload/MANIFEST.tsv  (file/size/sha256)
##   data-raw/zenodo_grouped_payload/gene_group_index.tsv
##   inst/extdata/gene_group_index.tsv             (shipped in package)

suppressMessages(devtools::load_all(".", quiet = TRUE))

OUT <- "data-raw/zenodo_grouped_payload"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## --- Load the real bundles ---------------------------------------
tp53  <- readRDS("tests/testthat/fixtures/TP53.rds")
wdr31 <- readRDS("tests/testthat/fixtures/WDR31.rds")
demo1 <- readRDS("inst/extdata/DEMO1.rds")

## --- Cheap synthetic bundles so we have >1 populated group --------
## empty_gene_data() leaves the ensembl id fields NA, which the strict
## validator rejects (no NA in required columns). Fill them with clearly
## synthetic placeholders so each bundle is spec-valid.
.make_synth <- function(gene, uniprot, len) {
    b <- empty_gene_data(gene, uniprot_id = uniprot, protein_length = len)
    b$meta$ensembl_gene_id       <- paste0("ENSG", gene)
    b$meta$ensembl_transcript_id <- paste0("ENST", gene)
    b$meta$source_versions       <- "SYNTHETIC — prototype only"
    b
}
synth1 <- .make_synth("SYNTH1", "P90001", 50L)
synth2 <- .make_synth("SYNTH2", "P90002", 60L)

## --- Pack groups: a group file is a NAMED LIST of gene bundles -----
group_01 <- list(TP53 = tp53, WDR31 = wdr31)
group_02 <- list(DEMO1 = demo1, SYNTH1 = synth1, SYNTH2 = synth2)

saveRDS(group_01, file.path(OUT, "group_01.rds"))
saveRDS(group_02, file.path(OUT, "group_02.rds"))
cat(sprintf("Wrote group_01 {%s} and group_02 {%s}\n",
            paste(names(group_01), collapse = ", "),
            paste(names(group_02), collapse = ", ")))

## --- Sanity: every packed bundle must validate -------------------
for (grp in list(group_01, group_02)) {
    for (g in names(grp)) {
        v <- validate_gene_data(grp[[g]])
        if (!isTRUE(v$valid)) {
            stop(sprintf("Bundle %s failed validation: %s",
                         g, paste(v$issues, collapse = "; ")))
        }
    }
}
cat("All packed bundles validate.\n")

## --- gene -> group index -----------------------------------------
index <- data.frame(
    gene  = c("TP53", "WDR31", "DEMO1", "SYNTH1", "SYNTH2"),
    group = c("group_01", "group_01", "group_02", "group_02", "group_02"),
    stringsAsFactors = FALSE
)
## Ships in the package (routing table, versioned with the package):
dir.create("inst/extdata", showWarnings = FALSE, recursive = TRUE)
write.table(index, "inst/extdata/gene_group_index.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
## Also drop a copy in the payload so the deposit is self-describing:
write.table(index, file.path(OUT, "gene_group_index.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("Wrote gene_group_index.tsv (5 genes -> 2 groups).\n")

## --- MANIFEST over the GROUP files (file/size/sha256) -------------
source("data-raw/build_scripts/build_manifest.R")
build_manifest(payload_dir = OUT)   # writes MANIFEST.tsv + README + LICENSES

cat("\nDone. Prototype grouped payload in", OUT, "\n")
