# Get gnomAD per-variant allele frequencies for a gene

Returns the \`gnomad\` slice with one row per coding variant in gnomAD
v4.1 (or whatever version was current at deposit build time; see
\`attr(., "source_versions")\`).

## Usage

``` r
get_gnomad(gene, force_refresh = FALSE)
```

## Arguments

- gene:

  HGNC gene symbol.

- force_refresh:

  Redownload even if cached.

## Value

A \`data.frame\` per \`DATA_FORMAT_SPEC.md\`. Returns \`NULL\` if
download failed.

## Examples

``` r
Sys.setenv(MSAVARIANT_CACHE = tempfile("msaVariant_cache_"))
import_local_bundle(
    system.file("extdata", "DEMO1.rds", package = "msaVariant"),
    gene = "DEMO1"
)
#> Cache directory does not exist yet (no annotations have been downloaded).
#> Imported DEMO1 bundle -> /tmp/Rtmp9JzdHU/msaVariant_cache_1a7351c13bd9/0.1.0/DEMO1.rds
head(get_gnomad("DEMO1"))
#>   pos aa_ref aa_alt aa_change consequence  af_joint ac_joint an_joint filter
#> 1   1      M      A       M1A    missense 0.0005818       12   152312   PASS
#> 2   5      Y      A       Y5A  synonymous 0.0009313       26   152312   PASS
#> 3   6      I      A       I6A    missense 0.0036489        1   152312   PASS
#> 4   8      K      A       K8A  synonymous 0.0020599        2   152312   PASS
#> 5   9      Q      A       Q9A    missense 0.0020708       24   152312   PASS
#> 6  11      Q      A      Q11A    missense 0.0024021       34   152312   PASS
```
