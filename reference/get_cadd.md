# Get CADD per-substitution PHRED scores for a gene

Per-substitution CADD PHRED-scaled deleteriousness scores. Higher = more
deleterious; PHRED 20 ~ top 1

## Usage

``` r
get_cadd(gene, force_refresh = FALSE)
```

## Arguments

- gene:

  HGNC gene symbol.

- force_refresh:

  Redownload even if cached.

## Value

A \`data.frame\` per spec; returns \`NULL\` on failure.

## Details

Note: CADD is licensed for \*\*non-commercial use only\*\*.

## Examples

``` r
Sys.setenv(MSAVARIANT_CACHE = tempfile("msaVariant_cache_"))
import_local_bundle(
    system.file("extdata", "DEMO1.rds", package = "msaVariant"),
    gene = "DEMO1"
)
#> Cache directory does not exist yet (no annotations have been downloaded).
#> Imported DEMO1 bundle -> /tmp/Rtmp9JzdHU/msaVariant_cache_1a735a4a0b17/0.1.0/DEMO1.rds
head(get_cadd("DEMO1"))
#>   pos aa_ref aa_alt aa_change consequence cadd_raw cadd_phred
#> 1   1      M      A       M1A    missense    1.535      12.28
#> 2   2      K      A       K2A    missense    0.539       4.31
#> 3   3      T      A       T3A    missense    4.897      39.17
#> 4   4      A      G       A4G    missense    2.485      19.88
#> 5   5      Y      A       Y5A    missense    0.465       3.72
#> 6   6      I      A       I6A    missense    1.059       8.47
```
