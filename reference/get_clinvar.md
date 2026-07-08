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
