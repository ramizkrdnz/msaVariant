# Overlay CADD per-residue PHRED scores onto an MSA

CADD is licensed for non-commercial use only.

## Usage

``` r
geom_cadd(
  data = NULL,
  gene = NULL,
  msa,
  ref_name = NULL,
  summary = c("median", "mean", "max"),
  value_range = c(0, 40),
  y_offset = -8,
  track_height = 1.2
)
```

## Arguments

- data:

  Optional pre-built data.frame; see \`DATA_FORMAT_SPEC.md\` for the
  schema. If \`NULL\`, supply \`gene\`.

- gene:

  HGNC symbol; triggers fetch from Zenodo deposit.

- msa, ref_name:

  Reference MSA and sequence name.

- summary:

  How to aggregate the 19 alt scores at each position. One of \`"mean"\`
  (default), \`"max"\`, \`"median"\`.

- value_range:

  Length-2 PHRED range for the colour scale. Default \`c(0, 40)\` covers
  the practically informative range.

- y_offset, track_height:

  Track geometry.

## Value

A list of ggplot2 layers (or \`NULL\`).
