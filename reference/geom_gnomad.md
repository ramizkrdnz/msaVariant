# Overlay gnomAD allele frequencies onto an MSA

Renders a per-residue heatmap strip of gnomAD allele frequencies. For
residues with multiple variants, the maximum joint allele frequency
across variants is plotted.

## Usage

``` r
geom_gnomad(
  data = NULL,
  gene = NULL,
  msa,
  ref_name = NULL,
  y_offset = -2,
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

- y_offset, track_height:

  Track geometry.

## Value

A list of ggplot2 layers (or \`NULL\`).
