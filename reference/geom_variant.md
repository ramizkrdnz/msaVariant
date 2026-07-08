# Overlay variants onto a multiple sequence alignment

Adds vertical bands marking the alignment columns affected by one or
more clinical / population variants. Works with a tidy data frame; for
VCF input use \`read_vcf_variants()\` first.

## Usage

``` r
geom_variant(
  variants,
  msa,
  ref_name = NULL,
  colour_by = c("consequence", "pathogenicity"),
  alpha = 0.5,
  show_label = TRUE,
  label_y = NULL,
  ...
)
```

## Arguments

- variants:

  A data.frame with at minimum a column \`pos\` (integer, 1-based
  protein residue position in the reference). Optional columns: \*
  \`pos_end\` : integer, for indels/frameshifts; defaults to \`pos\`
  (single-residue variant) \* \`label\` : character, shown above the
  band \* \`consequence\` : factor, used to colour the band (e.g.
  "missense", "frameshift", "nonsense", "synonymous", "splice") \*
  \`pathogenicity\` : factor, alternative colour key (e.g. "Pathogenic",
  "Likely_pathogenic", "VUS", "Likely_benign", "Benign")

- msa, ref_name:

  The MSA used in the parent \`ggmsa()\` call, and the reference
  sequence within it. The geom maps variant positions to alignment
  columns via \`map_variant_to_msa()\`.

- colour_by:

  Either \`"consequence"\` (default) or \`"pathogenicity"\`, controlling
  which column drives the fill.

- alpha:

  Band transparency in \`\[0, 1\]\`. Default \`0.5\`.

- show_label:

  Logical; whether to draw labels above the band.

- label_y:

  Numeric vertical position of the variant label. Defaults to \`n_seq +
  3\` (just above the MSA, clear of the default \`geom_domain()\` track
  at \`n_seq + 1.7\`).

- ...:

  Passed to the underlying \`geom_rect()\` layer.

## Value

A ggplot2 layer (or list of layers) that can be \`+\`-ed onto a
\`ggmsa()\` plot.

## Examples

``` r
if (FALSE) { # \dontrun{
library(ggmsa)
fa <- system.file("extdata", "patl1_orthologs.fasta",
                  package = "msaVariant")
v  <- data.frame(pos = 518, pos_end = 577,
                 label = "p.K518fs",
                 consequence = "frameshift")
ggmsa(fa, start = 500, end = 580) +
  geom_variant(v, msa = fa, ref_name = "PATL1_HUMAN")
} # }
```
