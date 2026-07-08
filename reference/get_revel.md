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
