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
  \[import_local_bundle()\] / \[fetch_gene_data()\]).

- variant_pos:

  Integer residue position (canonical UniProt numbering).

- aa_change:

  Character variant identifier, e.g. `"R175H"` or `"p.R175H"`. Ref/alt
  are parsed from this string.

- pp3_min_predictors:

  Integer 1–3 (default 2). PP3 fires when at least this many of the
  three PP3 predictors pass their threshold. `1` = any predictor, `2` =
  majority (default), `3` = all three must agree.

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

  `variant_pos` falls within the start–end range of at least one
  `domains` row.

- PM2:

  This `aa_change` is absent from gnomAD, or its `af_joint` is below
  0.0001.

- PM5:

  ClinVar has a record at the same position with a *different* alt
  residue and significance `"Pathogenic"` or `"Likely_pathogenic"`.

- PP3:

  At this substitution, count how many of the three computational
  predictors pass their threshold — AlphaMissense
  `am_class == "likely_pathogenic"`, REVEL `revel_score > 0.7`, CADD
  `cadd_phred >= 25`. Triggered when the passing count is at least
  `pp3_min_predictors`.

Codes whose evaluation requires the alternate allele (PS1, PM2, PM5,
PP3) are only evaluated when `aa_change` parses to a full ref+pos+alt
substitution (e.g. `"R175H"`). PM1 depends only on position and is
always evaluated. This keeps the function safe when called with a
position-only label (it simply won't over-fire PM2).

## Examples

``` r
## Runnable with the shipped synthetic DEMO1 bundle.
demo <- readRDS(system.file("extdata", "DEMO1.rds", package = "msaVariant"))
compute_acmg_codes(demo, variant_pos = 21, aa_change = "R21H")
#> [1] "PS1" "PM1" "PM2" "PM5" "PP3"

## The PP3 threshold is tunable: 1 = any predictor, 3 = all three.
compute_acmg_codes(demo, 30, "I30A", pp3_min_predictors = 3)
#> [1] "PM1"
```
