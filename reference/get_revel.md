# Get REVEL per-substitution scores for a gene

Per-substitution REVEL ensemble scores in \`\[0, 1\]\`. Higher indicates
more likely pathogenic.

## Usage

``` r
get_revel(gene, force_refresh = FALSE)
```

## Arguments

- gene:

  HGNC gene symbol.

- force_refresh:

  Redownload even if cached.

## Value

A \`data.frame\` per spec; returns \`NULL\` on failure.

## Examples

``` r
Sys.setenv(MSAVARIANT_CACHE = tempfile("msaVariant_cache_"))
import_local_bundle(
  system.file("extdata", "DEMO1.rds", package = "msaVariant"),
  gene = "DEMO1"
)
#> Cache directory does not exist yet (no annotations have been downloaded).
#> Imported DEMO1 bundle -> /tmp/RtmpUwZ6Qa/msaVariant_cache_1a0c75c08201/0.1.0/DEMO1.rds
head(get_revel("DEMO1"))
#>   pos aa_ref aa_alt aa_change revel_score
#> 1   1      M      A       M1A      0.8868
#> 2   2      K      A       K2A      0.1363
#> 3   3      T      A       T3A      0.7853
#> 4   4      A      G       A4G      0.4533
#> 5   5      Y      A       Y5A      0.1357
#> 6   6      I      A       I6A      0.8852
```
