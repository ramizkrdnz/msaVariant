# Compute triggered ACMG evidence codes for one variant

Evaluates a fixed subset of ACMG/AMP evidence codes for a single
missense substitution using only the annotation tables carried in a gene
bundle. Returns just the codes that are triggered.

## Usage

``` r
compute_acmg_codes(bundle, variant_pos, aa_change, pp3_min_predictors = 2L)
```

## Arguments

- bundle:

  A gene bundle: a named list with elements `domains`, `clinvar`,
  `gnomad`, `alphamissense`, `revel`, `cadd` (as produced by
  [`import_local_bundle`](https://ramizkrdnz.github.io/msaVariant/reference/import_local_bundle.md)
  /
  [`fetch_gene_data`](https://ramizkrdnz.github.io/msaVariant/reference/fetch_gene_data.md)).

- variant_pos:

  Integer residue position (canonical UniProt numbering).

- aa_change:

  Character variant identifier, e.g. `"R175H"` or `"p.R175H"`. Ref/alt
  are parsed from this string.

- pp3_min_predictors:

  Integer 1-3 (default 2). PP3 fires when at least this many of the
  three PP3 predictors (AlphaMissense class == "likely_pathogenic",
  REVEL \> 0.7, CADD phred \>= 25) pass their threshold.

## Value

A character vector of the triggered codes, in canonical
evidence-strength order (`PS1`, `PM1`, `PM2`, `PM5`, `PP3`);
`character(0)` if none fire.

## Details

Rules (all computed from the bundle; no external lookups):

- PS1:

  ClinVar contains the exact same `aa_change` with significance
  `"Pathogenic"`.

- PM1:

  `variant_pos` falls within the start-end range of at least one
  `domains` row.

- PM2:

  This `aa_change` is absent from gnomAD, or its `af_joint` is below
  0.0001.

- PM5:

  ClinVar has a record at the same position with a different alt residue
  and significance `"Pathogenic"` or `"Likely_pathogenic"`.

- PP3:

  Count how many of the three predictors pass — AlphaMissense
  `am_class == "likely_pathogenic"`, REVEL `revel_score > 0.7`, CADD
  `cadd_phred >= 25`; triggered when the count is at least
  `pp3_min_predictors`.

## Examples

``` r
if (FALSE) { # \dontrun{
b <- fetch_gene_data("TP53")
compute_acmg_codes(b, 175, "R175H")
} # }
```
