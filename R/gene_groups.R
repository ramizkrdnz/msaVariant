## -----------------------------------------------------------------
## msaVariant: grouped-bundle routing (Path 2 architecture)
## -----------------------------------------------------------------
##
## Instead of one .rds per gene, the data deposit packs many gene
## bundles into a smaller number of "group" .rds files (a group file
## deserialises to a NAMED LIST of gene bundles). This keeps the number
## of Zenodo files manageable for ~19,500 genes.
##
## A gene->group index (inst/extdata/gene_group_index.tsv, columns
## `gene` and `group`) ships in the package and routes a gene symbol to
## the group file that carries it. fetch_gene_data() downloads/caches
## the whole group once and extracts the requested gene's bundle, so a
## second gene in the same group needs no further download.
##
## Everything here degrades gracefully: if grouping is disabled or the
## gene is absent from the index, fetch_gene_data() falls back to the
## original per-gene mode.

## Grouped mode is on by default. Set MSAVARIANT_GROUPED to a false-y
## value ("0", "false", "no", "off") to force the legacy per-gene mode.
.grouped_enabled <- function() {
    v <- tolower(Sys.getenv("MSAVARIANT_GROUPED", "1"))
    !(v %in% c("0", "false", "no", "off", ""))
}

## Path to the gene->group index. Defaults to the copy shipped in the
## package; MSAVARIANT_GROUP_INDEX overrides it (used for testing a
## staged index before it is baked into the package).
.gene_group_index_path <- function() {
    ov <- Sys.getenv("MSAVARIANT_GROUP_INDEX", "")
    if (nzchar(ov)) {
        return(ov)
    }
    system.file("extdata", "gene_group_index.tsv", package = "msaVariant")
}

## Read the index as a data.frame with `gene` and `group`, or NULL if
## it is absent / unreadable / malformed.
.gene_group_index <- function() {
    p <- .gene_group_index_path()
    if (!nzchar(p) || !file.exists(p)) {
        return(NULL)
    }
    idx <- tryCatch(
        utils::read.delim(p,
            sep = "\t", stringsAsFactors = FALSE,
            colClasses = "character"
        ),
        error = function(e) NULL
    )
    if (is.null(idx) || !all(c("gene", "group") %in% names(idx))) {
        return(NULL)
    }
    idx
}

## Resolve a gene symbol to its group file basename, or NA_character_
## when there is no index / no entry for the gene.
.gene_group <- function(gene, index = .gene_group_index()) {
    if (is.null(index)) {
        return(NA_character_)
    }
    row <- index[index$gene == gene, , drop = FALSE]
    if (nrow(row) < 1L) {
        return(NA_character_)
    }
    grp <- row$group[1]
    if (is.null(grp) || is.na(grp) || !nzchar(grp)) {
        return(NA_character_)
    }
    grp
}

## Cache path for a group file: <cache>/<data version>/<group>.rds,
## alongside any per-gene bundles.
.group_cache_path <- function(group) {
    d <- file.path(.cache_dir(), MSAVARIANT_DATA_VERSION)
    if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
    file.path(d, paste0(group, ".rds"))
}
