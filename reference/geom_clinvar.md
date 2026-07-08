# Overlay ClinVar pathogenicity calls onto an MSA

Renders a discrete tick row showing ClinVar calls at each affected
residue. Where multiple ClinVar entries exist at one residue, the most
pathogenic call is shown.

## Usage

``` r
geom_clinvar(
  data = NULL,
  gene = NULL,
  msa,
  ref_name = NULL,
  significance = c("Pathogenic", "Likely_pathogenic", "VUS", "Likely_benign", "Benign"),
  y_offset = -3.5,
  track_height = 1
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

- significance:

  Character vector of significance levels to show (default: all five
  ACMG tiers).

- y_offset, track_height:

  Track geometry.

## Value

A list of ggplot2 layers (or \`NULL\`).
