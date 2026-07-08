# Get gene metadata

Returns the \`meta\` slice — one row with UniProt ID, protein length,
Ensembl IDs, build date, and source-version stamps.

## Usage

``` r
get_gene_meta(gene, force_refresh = FALSE)
```

## Arguments

- gene:

  HGNC gene symbol.

- force_refresh:

  Redownload even if cached.

## Value

A 1-row \`data.frame\`; returns \`NULL\` on failure.
