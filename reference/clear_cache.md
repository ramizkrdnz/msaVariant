# Clear the msaVariant download cache

Removes locally-cached annotation files. Use this if you suspect a
cached file is stale or corrupt, or if you want to free disk space.

## Usage

``` r
clear_cache(gene = NULL)
```

## Arguments

- gene:

  If supplied, only that gene's cached file is removed. If \`NULL\`
  (default), the entire cache is removed.

## Value

Invisibly, the number of files deleted.

## Examples

``` r
Sys.setenv(MSAVARIANT_CACHE = tempfile("msaVariant_cache_"))
import_local_bundle(
  system.file("extdata", "DEMO1.rds", package = "msaVariant"),
  gene = "DEMO1"
)
#> Cache directory does not exist yet (no annotations have been downloaded).
#> Imported DEMO1 bundle -> /tmp/RtmpUwZ6Qa/msaVariant_cache_1a0c3824b5b9/0.1.0/DEMO1.rds
clear_cache("DEMO1")   # remove one gene
#> Cleared cached file for DEMO1.
clear_cache()          # remove everything
#> Cleared 0 cached file(s).
```
