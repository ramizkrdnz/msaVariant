# Get ClinVar variants for a gene

Returns the \`clinvar\` slice. See \`DATA_FORMAT_SPEC.md\` for column
definitions.

## Usage

``` r
get_clinvar(gene, force_refresh = FALSE)
```

## Arguments

- gene:

  HGNC gene symbol.

- force_refresh:

  Redownload even if cached.

## Value

A \`data.frame\` with \`pos\`, \`aa_ref\`, \`aa_alt\`, \`aa_change\`,
\`significance\` (factor), \`review_status\` (factor), \`clinvar_id\`,
and optionally \`condition\`, \`last_evaluated\`. Returns \`NULL\` if
download failed.

## Examples

``` r
Sys.setenv(MSAVARIANT_CACHE = tempfile("msaVariant_cache_"))
import_local_bundle(
  system.file("extdata", "DEMO1.rds", package = "msaVariant"),
  gene = "DEMO1"
)
#> Cache directory does not exist yet (no annotations have been downloaded).
#> Imported DEMO1 bundle -> /tmp/RtmpUwZ6Qa/msaVariant_cache_1a0c41b5f75b/0.1.0/DEMO1.rds
head(get_clinvar("DEMO1"))
#>   pos aa_ref aa_alt aa_change      significance review_status  clinvar_id
#> 1   1      M      A       M1A Likely_pathogenic        2_star DEMO0000001
#> 2   2      K      A       K2A             Other        3_star DEMO0000002
#> 3   3      T      A       T3A     Likely_benign        1_star DEMO0000003
#> 4   4      A      G       A4G       Conflicting        3_star DEMO0000004
#> 5   5      Y      A       Y5A            Benign        1_star DEMO0000005
#> 6   7      A      G       A7G            Benign        1_star DEMO0000006
```
