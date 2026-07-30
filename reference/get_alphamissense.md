# Get AlphaMissense per-substitution scores for a gene

Returns the full per-substitution table (every missense substitution at
every residue, ~19 rows per residue).

## Usage

``` r
get_alphamissense(gene, force_refresh = FALSE)
```

## Arguments

- gene:

  HGNC gene symbol.

- force_refresh:

  Redownload even if cached.

## Value

A \`data.frame\` with \`pos\`, \`aa_ref\`, \`aa_alt\`, \`aa_change\`,
\`am_score\`, \`am_class\` (factor: likely_benign / ambiguous /
likely_pathogenic). Returns \`NULL\` if download failed.

## Details

Note: AlphaMissense is licensed CC-BY-NC-SA 4.0. Commercial use is
restricted.

## Examples

``` r
Sys.setenv(MSAVARIANT_CACHE = tempfile("msaVariant_cache_"))
import_local_bundle(
    system.file("extdata", "DEMO1.rds", package = "msaVariant"),
    gene = "DEMO1"
)
#> Cache directory does not exist yet (no annotations have been downloaded).
#> Imported DEMO1 bundle -> /tmp/Rtmp9JzdHU/msaVariant_cache_1a7366a51d6f/0.1.0/DEMO1.rds
head(get_alphamissense("DEMO1"))
#>   pos aa_ref aa_alt aa_change am_score          am_class
#> 1   1      M      A       M1A   0.2457     likely_benign
#> 2   2      K      A       K2A   0.3511         ambiguous
#> 3   3      T      A       T3A   0.1590     likely_benign
#> 4   4      A      G       A4G   0.3041     likely_benign
#> 5   5      Y      A       Y5A   0.0175     likely_benign
#> 6   6      I      A       I6A   0.9966 likely_pathogenic
```
