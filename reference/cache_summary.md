# Summarize what's currently cached

Lists all locally-cached per-gene annotation files, with their size and
date of caching. Useful for inspecting disk use after querying many
genes.

## Usage

``` r
cache_summary()
```

## Value

A \`data.frame\` with columns \`gene\`, \`size_kb\`, \`cached_on\`,
sorted by size descending. Returns an empty data.frame if the cache is
empty.

## Examples

``` r
Sys.setenv(MSAVARIANT_CACHE = tempfile("msaVariant_cache_"))
import_local_bundle(
    system.file("extdata", "DEMO1.rds", package = "msaVariant"),
    gene = "DEMO1"
)
#> Cache directory does not exist yet (no annotations have been downloaded).
#> Imported DEMO1 bundle -> /tmp/Rtmp9JzdHU/msaVariant_cache_1a73121b0561/0.1.0/DEMO1.rds
cache_summary()
#>    gene size_kb  cached_on
#> 1 DEMO1     3.1 2026-07-30
```
