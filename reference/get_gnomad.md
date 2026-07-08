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
