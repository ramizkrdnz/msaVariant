## Tests for the grouped-bundle fetch architecture (Path 2).
## No network: a staging dir holds the group .rds files and
## MSAVARIANT_LOCAL_SOURCE makes the "download" a local copy, so we can
## prove download-once + cache-hit behaviour deterministically.

# Build a staging deposit (group_01 = {TP53, WDR31}) plus a matching
# gene->group index, and point the package's env vars at them. Returns
# the staging + cache dirs. All dirs are withr-scoped to the caller.
.setup_grouped <- function(env = parent.frame()) {
    stage <- withr::local_tempdir(.local_envir = env)
    cache <- withr::local_tempdir(.local_envir = env)
    tp53 <- readRDS(test_path("fixtures", "TP53.rds"))
    wdr31 <- readRDS(test_path("fixtures", "WDR31.rds"))
    saveRDS(list(TP53 = tp53, WDR31 = wdr31), file.path(stage, "group_01.rds"))
    idx <- data.frame(
        gene = c("TP53", "WDR31"),
        group = c("group_01", "group_01"),
        stringsAsFactors = FALSE
    )
    idxfile <- file.path(stage, "gene_group_index.tsv")
    write.table(idx, idxfile, sep = "\t", row.names = FALSE, quote = FALSE)
    # A MANIFEST.tsv over the group file, so the fetcher's manifest
    # auto-download (served from the staging dir) has something to fetch.
    grp <- file.path(stage, "group_01.rds")
    write.table(
        data.frame(
            file = "group_01.rds",
            size = file.info(grp)$size,
            sha256 = .sha256_file(grp),
            stringsAsFactors = FALSE
        ),
        file.path(stage, "MANIFEST.tsv"),
        sep = "\t", row.names = FALSE, quote = FALSE
    )
    withr::local_envvar(
        MSAVARIANT_CACHE = cache,
        MSAVARIANT_LOCAL_SOURCE = stage,
        MSAVARIANT_GROUP_INDEX = idxfile,
        MSAVARIANT_GROUPED = "1",
        .local_envir = env
    )
    list(stage = stage, cache = cache, idxfile = idxfile)
}

test_that("index routes a gene to its group", {
    s <- .setup_grouped()
    expect_equal(.gene_group("TP53"), "group_01")
    expect_equal(.gene_group("WDR31"), "group_01")
    expect_true(is.na(.gene_group("NOSUCHGENE")))
})

test_that("grouped fetch downloads the group once and extracts the gene", {
    s <- .setup_grouped()
    b <- fetch_gene_data("TP53", quiet = TRUE)
    expect_type(b, "list")
    expect_length(b, 7L)
    expect_equal(b$meta$gene, "TP53")
    expect_true(validate_gene_data(b)$valid)

    # The GROUP file is cached; no per-gene file is written.
    grp <- file.path(s$cache, "0.1.0", "group_01.rds")
    expect_true(file.exists(grp))
    expect_false(file.exists(file.path(s$cache, "0.1.0", "TP53.rds")))
})

test_that("second gene in the same group is served from cache (no re-download)", {
    s <- .setup_grouped()
    b1 <- fetch_gene_data("TP53", quiet = TRUE)
    expect_equal(b1$meta$gene, "TP53")

    # Delete the source: a genuine re-download would now fail. A cache hit
    # on the already-downloaded group still succeeds.
    file.remove(file.path(s$stage, "group_01.rds"))
    b2 <- fetch_gene_data("WDR31", quiet = TRUE)
    expect_equal(b2$meta$gene, "WDR31")
    expect_true(validate_gene_data(b2)$valid)
})

test_that("MSAVARIANT_GROUPED=0 forces legacy per-gene mode", {
    s <- .setup_grouped()
    withr::local_envvar(MSAVARIANT_GROUPED = "0")
    # Flag off: the group index is bypassed, so TP53 is fetched per-gene.
    # Stage a per-gene TP53.rds so the per-gene path can serve it.
    file.copy(
        test_path("fixtures", "TP53.rds"),
        file.path(s$stage, "TP53.rds")
    )
    b <- fetch_gene_data("TP53", quiet = TRUE)
    expect_equal(b$meta$gene, "TP53")
    # Per-gene mode writes a per-gene cache file, not a group file.
    expect_true(file.exists(file.path(s$cache, "0.1.0", "TP53.rds")))
    expect_false(file.exists(file.path(s$cache, "0.1.0", "group_01.rds")))
})

test_that("an imported per-gene bundle takes precedence over grouped mode", {
    s <- .setup_grouped()
    # Seed a per-gene bundle directly in the cache.
    dir.create(file.path(s$cache, "0.1.0"), recursive = TRUE, showWarnings = FALSE)
    demo <- readRDS(system.file("extdata", "DEMO1.rds", package = "msaVariant"))
    saveRDS(demo, file.path(s$cache, "0.1.0", "TP53.rds")) # deliberately mislabelled
    b <- fetch_gene_data("TP53", quiet = TRUE)
    # The per-gene cache file wins, so we get DEMO1's content, not group_01's.
    expect_equal(b$meta$gene, "DEMO1")
})

test_that("a corrupt cached group is rejected by checksum verification", {
    s <- .setup_grouped()
    # Write a MANIFEST listing the true sha256 of the group file.
    grp_src <- file.path(s$stage, "group_01.rds")
    manifest <- data.frame(
        file = "group_01.rds",
        size = file.info(grp_src)$size,
        sha256 = .sha256_file(grp_src),
        stringsAsFactors = FALSE
    )
    verdir <- file.path(s$cache, "0.1.0")
    dir.create(verdir, recursive = TRUE, showWarnings = FALSE)
    write.table(manifest, file.path(verdir, "MANIFEST.tsv"),
        sep = "\t", row.names = FALSE, quote = FALSE
    )
    # First fetch downloads + caches the group (checksum matches).
    expect_equal(fetch_gene_data("TP53", quiet = TRUE)$meta$gene, "TP53")

    # Corrupt the cached group, then point the source at a bad file too so
    # a re-download can't silently repair it: the corrupt cache must be
    # detected and rejected.
    grp_cache <- file.path(verdir, "group_01.rds")
    writeBin(as.raw(rep(0L, 200)), grp_cache)
    file.remove(grp_src) # no valid source to re-download from
    w <- capture_warnings(out <- fetch_gene_data("TP53", quiet = TRUE))
    expect_true(any(grepl("checksum", w)))
    expect_null(out)
})

test_that("fetch auto-downloads MANIFEST.tsv into the cache once", {
    s <- .setup_grouped() # staging dir includes a MANIFEST.tsv
    man_cache <- file.path(s$cache, "0.1.0", "MANIFEST.tsv")
    expect_false(file.exists(man_cache)) # not there before any fetch

    b <- fetch_gene_data("TP53", quiet = TRUE)
    expect_equal(b$meta$gene, "TP53")
    # The manifest was pulled automatically alongside the group file.
    expect_true(file.exists(man_cache))
    man <- utils::read.delim(man_cache, sep = "\t", stringsAsFactors = FALSE)
    expect_true("group_01.rds" %in% man$file)
})

test_that("checksum verification engages automatically via auto-download", {
    # No manifest is placed manually: the fetcher must pull it itself and
    # then reject a corrupt cached group on the next call.
    s <- .setup_grouped()
    expect_equal(fetch_gene_data("TP53", quiet = TRUE)$meta$gene, "TP53")
    verdir <- file.path(s$cache, "0.1.0")
    expect_true(file.exists(file.path(verdir, "MANIFEST.tsv"))) # auto-pulled

    # Corrupt the cached group; remove the source so it can't be repaired.
    writeBin(as.raw(rep(0L, 200)), file.path(verdir, "group_01.rds"))
    file.remove(file.path(s$stage, "group_01.rds"))
    w <- capture_warnings(out <- fetch_gene_data("TP53", quiet = TRUE))
    expect_true(any(grepl("checksum", w)))
    expect_null(out)
})

test_that("missing manifest degrades gracefully (fetch still succeeds)", {
    # Remove the staged manifest: auto-download finds nothing, so
    # verification stays dormant and the fetch still returns the bundle.
    s <- .setup_grouped()
    file.remove(file.path(s$stage, "MANIFEST.tsv"))
    b <- fetch_gene_data("TP53", quiet = TRUE)
    expect_equal(b$meta$gene, "TP53")
    expect_true(validate_gene_data(b)$valid)
    # No manifest landed in the cache; behaviour is as before.
    expect_false(file.exists(file.path(s$cache, "0.1.0", "MANIFEST.tsv")))
})
