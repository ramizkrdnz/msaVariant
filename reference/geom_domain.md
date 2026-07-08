# Overlay protein domain annotation on an MSA

Draws a horizontal track showing protein domain ranges with labels,
mapped to MSA columns via the reference sequence.

## Usage

``` r
geom_domain(
  domains = NULL,
  gene = NULL,
  msa,
  ref_name = NULL,
  y_offset = NULL,
  track_height = 1,
  fill = "#999999"
)
```

## Arguments

- domains:

  data.frame with columns \`start\`, \`end\`, \`name\` (residue
  positions, 1-based, inclusive). Optional column \`source\` (e.g.
  "InterPro", "Pfam", "SMART") is preserved. If \`NULL\`, supply
  \`gene\` instead.

- gene:

  HGNC gene symbol; if supplied, the package will fetch domain
  annotations from the Zenodo deposit.

- msa, ref_name:

  Reference MSA and sequence.

- y_offset, track_height:

  Track geometry.

- fill:

  Default fill colour for domain boxes if no \`domains\$colour\` column
  is provided.

## Value

A list of ggplot2 layers.
