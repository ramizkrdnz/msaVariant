# Overlay REVEL per-residue scores onto an MSA

Overlay REVEL per-residue scores onto an MSA

## Usage

``` r
geom_revel(
  data = NULL,
  gene = NULL,
  msa,
  ref_name = NULL,
  summary = c("mean", "max", "median"),
  y_offset = -6.5,
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

- y_offset, track_height:

  Track geometry.

## Value

A list of ggplot2 layers (or \`NULL\`).
