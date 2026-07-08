# Overlay an arbitrary per-residue annotation track onto an MSA

Renders a strip below (or above) the MSA showing any per-residue scalar
or categorical value the user supplies. Handles both continuous data
(rendered as a colour-graded heatmap strip) and discrete data (rendered
as a coloured tick row).

## Usage

``` r
geom_track(
  data,
  msa,
  ref_name = NULL,
  value,
  type = c("continuous", "discrete"),
  name = "",
  palette = NULL,
  value_range = NULL,
  y_offset = -2,
  track_height = 1.2
)
```

## Arguments

- data:

  A data.frame with at least: \* \`pos\` — integer residue position
  (1-based in the reference) \* a value column whose name you pass via
  \`value\`

- msa, ref_name:

  MSA and reference sequence name (same as in the parent \`ggmsa()\`
  call).

- value:

  Name of the column in \`data\` holding the per-residue value to plot.

- type:

  Either \`"continuous"\` (heatmap strip, default) or \`"discrete"\`
  (categorical tick row).

- name:

  Label printed to the right of the track.

- palette:

  For \`type = "continuous"\`, a vector of colours passed to
  \`colorRampPalette()\`. For \`type = "discrete"\`, a named vector
  mapping levels of \`value\` to colours.

- value_range:

  For continuous tracks, a length-2 numeric vector giving \`c(min,
  max)\` for the colour scale. Defaults to the observed range.

- y_offset, track_height:

  Vertical placement and height in MSA-row units.

## Value

A list of ggplot2 layers.

## Details

This is the generic, bring-your-own-data layer for annotation tracks.
Use it for gnomAD allele frequencies, ClinVar significance,
AlphaMissense / REVEL / CADD scores, post-translational-modification
sites, ChIP-seq signal, or anything else you can map to a residue.

## Examples

``` r
if (FALSE) { # \dontrun{
# Continuous: gnomAD allele frequencies from your own query
my_af <- data.frame(pos = 510:530, af = runif(21, 0, 1e-3))
geom_track(my_af, msa = fa, value = "af",
           name = "gnomAD AF", type = "continuous")

# Discrete: ClinVar calls from your own export
my_cv <- data.frame(pos = c(498, 518),
                    sig = c("Benign", "Pathogenic"))
geom_track(my_cv, msa = fa, value = "sig", type = "discrete",
           name = "ClinVar",
           palette = c(Benign = "#4575B4", Pathogenic = "#D7301F"))
} # }
```
