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
