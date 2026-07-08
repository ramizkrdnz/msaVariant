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
