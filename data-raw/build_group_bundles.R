## =================================================================
## build_group_bundles.R -- size-balanced grouped-bundle packer
## =================================================================
##
## Generalizes the 5-gene build_group_prototype.R into the real packer.
## Consumes a directory of per-gene <GENE>.rds bundles (the output of
## merge_into_gene_bundles) and packs them into group files for the
## Path-2 grouped deposit, then writes the gene->group index and a
## MANIFEST over the group files.
##
## Packing: alphabetical greedy with a DUAL cap -- a new group starts
## when adding the next gene would exceed EITHER `max_genes` (~200) OR
## `max_bytes` (~40 MB). Deterministic (genes sorted by symbol), so the
## same input always yields the same grouping. A single gene larger than
## `max_bytes` (e.g. TTN) lands in its own group rather than looping.
##
## Outputs into `out_dir`:
##   group_001.rds ... group_NNN.rds   (each a named list of gene bundles)
##   gene_group_index.tsv              (gene -> group, self-describing copy)
##   MANIFEST.tsv, README.md, LICENSES/ (via build_manifest)
## and ships the index at `index_dest` (default inst/extdata/), which is
## the routing table the package reads at runtime.
##
## USAGE (from package root):
##   suppressMessages(devtools::load_all("."))          # for validate_gene_data
##   source("data-raw/build_scripts/build_manifest.R")
##   source("data-raw/build_group_bundles.R")
##   build_group_bundles(
##     gene_dir = "zenodo_payload",                     # merge output
##     out_dir  = "zenodo_grouped_payload"
##   )

build_group_bundles <- function(gene_dir,
                                 out_dir,
                                 index_dest = "inst/extdata/gene_group_index.tsv",
                                 max_genes = 200L,
                                 max_bytes = 40L * 1024L^2,
                                 validate = TRUE,
                                 make_manifest = TRUE) {
  stopifnot(dir.exists(gene_dir))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  files <- list.files(gene_dir, pattern = "\\.rds$", full.names = FALSE)
  # Exclude any non-gene artefacts that might share the dir.
  files <- files[!files %in% c("MANIFEST.tsv")]
  if (length(files) == 0L) stop("No <gene>.rds files found in ", gene_dir)
  genes <- sub("\\.rds$", "", files)
  ord   <- order(genes)                    # deterministic, alphabetical
  genes <- genes[ord]; files <- files[ord]
  sizes <- file.info(file.path(gene_dir, files))$size

  ## ---- Greedy dual-cap packing into groups -----------------------
  group_of <- integer(length(genes))       # group index per gene
  g <- 1L; cur_n <- 0L; cur_bytes <- 0
  for (i in seq_along(genes)) {
    would_exceed <- (cur_n >= max_genes) || (cur_bytes + sizes[i] > max_bytes)
    if (would_exceed && cur_n > 0L) {       # close current group, start new
      g <- g + 1L; cur_n <- 0L; cur_bytes <- 0
    }
    group_of[i] <- g
    cur_n <- cur_n + 1L; cur_bytes <- cur_bytes + sizes[i]
  }
  n_groups <- g
  width <- max(3L, nchar(as.character(n_groups)))
  group_name <- sprintf(paste0("group_%0", width, "d"), group_of)

  ## ---- Write each group file + collect index rows ----------------
  index <- data.frame(gene = character(0), group = character(0),
                      stringsAsFactors = FALSE)
  excluded <- character(0)
  uniq_groups <- unique(group_name)
  for (gn in uniq_groups) {
    idx <- which(group_name == gn)
    bundles <- list()
    for (j in idx) {
      b <- readRDS(file.path(gene_dir, files[j]))
      if (isTRUE(validate)) {
        v <- validate_gene_data(b)
        if (!isTRUE(v$valid)) {
          warning(sprintf("Excluding %s from %s (failed validation): %s",
                          genes[j], gn, paste(v$issues, collapse = "; ")))
          excluded <- c(excluded, genes[j])
          next
        }
      }
      bundles[[genes[j]]] <- b
    }
    if (length(bundles) == 0L) next
    save_atomic(bundles, file.path(out_dir, paste0(gn, ".rds")))
    index <- rbind(index, data.frame(
      gene = names(bundles), group = gn, stringsAsFactors = FALSE))
  }
  index <- index[order(index$gene), , drop = FALSE]
  rownames(index) <- NULL

  ## ---- Write the gene->group index (shipped + self-describing) ---
  dir.create(dirname(index_dest), showWarnings = FALSE, recursive = TRUE)
  write.table(index, index_dest, sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(index, file.path(out_dir, "gene_group_index.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  ## ---- MANIFEST over the group files -----------------------------
  if (isTRUE(make_manifest)) {
    if (!exists("build_manifest", mode = "function")) {
      stop("build_manifest() not found; source data-raw/build_scripts/build_manifest.R first.")
    }
    build_manifest(payload_dir = out_dir)
  }

  ## ---- Summary ---------------------------------------------------
  grp_sizes <- table(index$group)
  message(sprintf(
    "Packed %d genes into %d groups (%d excluded). Group sizes: %d..%d genes.",
    nrow(index), length(unique(index$group)), length(excluded),
    min(grp_sizes), max(grp_sizes)))
  message(sprintf("  Index: %s (+ copy in %s)", index_dest, out_dir))
  invisible(list(index = index, n_groups = length(unique(index$group)),
                 excluded = excluded))
}

## Atomic write helper (shared shape with the other builders): write to a
## temp file in the same dir, then rename, so a crash never leaves a
## truncated group file.
save_atomic <- function(obj, dest) {
  tmp <- tempfile(tmpdir = dirname(dest), fileext = ".rds.part")
  saveRDS(obj, tmp, compress = "xz")
  if (!file.rename(tmp, dest)) {
    file.copy(tmp, dest, overwrite = TRUE); unlink(tmp)
  }
}
