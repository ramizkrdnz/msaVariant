## =================================================================
## build_manifest.R
## =================================================================

build_manifest <- function(payload_dir) {
  out_dir <- payload_dir
  files <- list.files(out_dir, pattern = "\\.rds$", full.names = TRUE)
  rel   <- basename(files)
  sizes <- file.info(files)$size

  message("  Computing SHA256 checksums ...")
  if (requireNamespace("digest", quietly = TRUE)) {
    shas <- vapply(files, function(p)
      digest::digest(file = p, algo = "sha256"),
      character(1L))
  } else {
    shas <- vapply(files, function(p) {
      out <- system2("sha256sum", p, stdout = TRUE)
      sub(" .*", "", out)
    }, character(1L))
  }

  manifest <- data.frame(
    file   = rel,
    size   = sizes,
    sha256 = unname(shas),
    stringsAsFactors = FALSE
  )
  write.table(manifest, file.path(out_dir, "MANIFEST.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  message(sprintf("  MANIFEST.tsv: %d entries", nrow(manifest)))

  ## ---- Licenses --------------------------------------------------
  lic_dir <- file.path(out_dir, "LICENSES")
  dir.create(lic_dir, showWarnings = FALSE)
  writeLines(c(
    "AlphaMissense scores (Cheng et al., Science 2023) are",
    "licensed under Creative Commons Attribution-NonCommercial-",
    "ShareAlike 4.0 International (CC-BY-NC-SA 4.0).",
    "",
    "See https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode"
  ), file.path(lic_dir, "AlphaMissense_LICENSE.txt"))

  writeLines(c(
    "CADD (Combined Annotation Dependent Depletion) is licensed",
    "for non-commercial use only by the University of Washington.",
    "",
    "Citation: Schubach M et al. NAR (2024) 'CADD v1.7'.",
    "doi:10.1093/nar/gkad989"
  ), file.path(lic_dir, "CADD_LICENSE.txt"))

  writeLines(c(
    "ClinVar data are made available by the National Center for",
    "Biotechnology Information (NCBI). NCBI does not specifically",
    "endorse any commercial product. ClinVar data are in the",
    "public domain.",
    "",
    "Citation: Landrum MJ et al. NAR (2020) 48:D835."
  ), file.path(lic_dir, "ClinVar_DISCLAIMER.txt"))

  ## ---- README ----------------------------------------------------
  readme <- c(
    "# msaVariant data deposit v0.1.0",
    "",
    sprintf("Build date: %s", Sys.Date()),
    "",
    "Per-gene annotation bundles for human protein-coding genes,",
    "consumed by the `msaVariant` R package",
    "(https://github.com/KaplanLab/msaVariant).",
    "",
    "## Format",
    "",
    "One .rds file per HGNC-approved gene symbol. Each file",
    "deserializes to an R list with elements: meta, domains,",
    "clinvar, gnomad, alphamissense, revel, cadd.",
    "See DATA_FORMAT_SPEC.md in the package repository for the",
    "complete schema.",
    "",
    "## Sources and snapshot dates",
    "",
    "| Source | Snapshot date | Citation |",
    "|---|---|---|",
    "| InterPro/Pfam | (build date) | Paysan-Lafosse T et al. NAR (2023) |",
    "| ClinVar       | (build date) | Landrum MJ et al. NAR (2020) |",
    "| gnomAD v4.1   | (build date) | Chen S et al. Nature (2024) 625:92 |",
    "| AlphaMissense | 2023-09      | Cheng J et al. Science (2023) 381 |",
    "| REVEL v1.3    | (build date) | Ioannidis NM et al. AJHG (2016) 99:877 |",
    "| CADD v1.7     | (build date) | Schubach M et al. NAR (2024) |",
    "",
    "## License",
    "",
    "This deposit aggregates data with mixed licensing. The most",
    "restrictive applicable license (CC-BY-NC-SA from AlphaMissense,",
    "plus the non-commercial restriction from CADD) propagates to",
    "the deposit as a whole.",
    "",
    "**Academic use: permitted.**",
    "**Commercial use: requires separate licenses from",
    "Google DeepMind (AlphaMissense) and the University of Washington",
    "(CADD).**",
    "",
    "## Reproducibility",
    "",
    "The build scripts that produced this deposit are in the package",
    "repository under data-raw/. To reproduce, install the package,",
    "then run `Rscript data-raw/build_gene_files.R`."
  )
  writeLines(readme, file.path(out_dir, "README.md"))
  message("  README.md, LICENSES/ written")
}
