# Per-column conservation score for an MSA

Per-column conservation score for an MSA

## Usage

``` r
conservation_score(
  msa,
  method = c("shannon", "js"),
  background = NULL,
  ignore_gaps = TRUE
)
```

## Arguments

- msa:

  MSA in any form accepted by \`build_msa_coord_map()\`.

- method:

  One of \`"shannon"\` (default) or \`"js"\` (Jensen-Shannon divergence
  against a background).

- background:

  Optional named vector of background residue frequencies. If \`NULL\`
  and \`method = "js"\`, a uniform background is used.

- ignore_gaps:

  If \`TRUE\` (default), gap characters are excluded from the per-column
  frequency calculation.

## Value

A data.frame with columns \`msa_col\` and \`score\`. For Shannon entropy
the score is in nats, rescaled so 1 = fully conserved and 0 = uniform
across the observed alphabet. For Jensen-Shannon, score is in \[0, 1\]
with 1 = fully conserved.
