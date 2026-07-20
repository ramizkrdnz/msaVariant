# Toggling Layers

A figure for a supervisor and a figure for a journal are rarely the same
figure.
[`plot_variant_overlay()`](https://ramizkrdnz.github.io/msaVariant/reference/plot_variant_overlay.md)
gives every annotation layer an on/off switch, and the layout re-flows
to whatever remains — no blank rows, no orphaned axis, no legend for a
track that is not there.

``` r

library(msaVariant)

Sys.setenv(MSAVARIANT_CACHE = tempfile("msaVariant_cache_"))
import_local_bundle(
  system.file("extdata", "DEMO1.rds", package = "msaVariant"),
  gene = "DEMO1"
)
demo_fasta <- system.file("extdata", "demo_aligned.fasta",
                          package = "msaVariant")
```

> DEMO1 is synthetic demonstration data, not real annotations.

## The switches

| Argument             | Layer                         |
|----------------------|-------------------------------|
| `show_acmg`          | ACMG evidence-code chip strip |
| `show_domain`        | Protein domain track          |
| `show_clinvar`       | ClinVar significance track    |
| `show_gnomad`        | gnomAD allele-frequency track |
| `show_cadd`          | CADD PHRED track              |
| `show_alphamissense` | AlphaMissense score track     |
| `show_revel`         | REVEL score track             |

All default to `TRUE`. The variant header and the MSA panel have no
switch — they are the figure.

## Everything on

The default: header, ACMG strip, alignment, all six tracks, and a legend
dashboard.

``` r

plot_variant_overlay(
  gene          = "DEMO1",
  aligned_fasta = demo_fasta,
  variant_pos   = 21,
  variant_label = "p.R21H"
)
```

![Full DEMO1 p.R21H figure with all annotation
layers](toggling-layers_files/figure-html/full-1.png)

## Everything off but the essentials

Alignment, variant, ACMG strip. Nothing else.

``` r

plot_variant_overlay(
  gene               = "DEMO1",
  aligned_fasta      = demo_fasta,
  variant_pos        = 21,
  variant_label      = "p.R21H",
  show_domain        = FALSE,
  show_clinvar       = FALSE,
  show_gnomad        = FALSE,
  show_cadd          = FALSE,
  show_alphamissense = FALSE,
  show_revel         = FALSE
)
```

![Minimal DEMO1 p.R21H figure: MSA, variant and ACMG strip
only](toggling-layers_files/figure-html/minimal-1.png)

Compare the two. The minimal figure is not the full figure with holes in
it — the panel heights are recomputed, the residue-position axis has
migrated from the bottom track up onto the MSA panel, and the legend
dashboard has disappeared entirely because no legend-carrying track
remains. That is the whole point of the toggles: a subset stays
publication-ready.

## In-between: the in-silico predictors only

A common case — you want the computational evidence and the alignment,
without the clinical and population context.

``` r

plot_variant_overlay(
  gene          = "DEMO1",
  aligned_fasta = demo_fasta,
  variant_pos   = 21,
  variant_label = "p.R21H",
  show_domain   = FALSE,
  show_clinvar  = FALSE,
  show_gnomad   = FALSE
)
```

![DEMO1 p.R21H figure with only the in-silico predictor
tracks](toggling-layers_files/figure-html/predictors-1.png)

Three tracks remain, three legends appear, and the axis sits on REVEL —
the lowest surviving track.

## Dropping the ACMG strip

`show_acmg = FALSE` removes the chip row. The codes are still computed
internally; they are simply not drawn.

``` r

plot_variant_overlay(
  gene          = "DEMO1",
  aligned_fasta = demo_fasta,
  variant_pos   = 21,
  variant_label = "p.R21H",
  show_acmg     = FALSE
)
```

![DEMO1 p.R21H figure with the ACMG strip
hidden](toggling-layers_files/figure-html/no-acmg-1.png)

## Subsetting the codes instead

If you want the strip but not all five chips, keep `show_acmg = TRUE`
and pass `acmg_codes` — the display filter. Only codes that are both
requested *and* triggered are drawn, so asking for a code that did not
fire simply shows nothing rather than an empty chip.

``` r

plot_variant_overlay(
  gene          = "DEMO1",
  aligned_fasta = demo_fasta,
  variant_pos   = 21,
  variant_label = "p.R21H",
  acmg_codes    = c("PS1", "PM5")
)
```

![DEMO1 p.R21H figure showing only the PS1 and PM5
chips](toggling-layers_files/figure-html/acmg-subset-1.png)

The distinction is worth holding onto:

- `show_acmg` — does the strip exist at all?
- `acmg_codes` — which of the triggered chips get drawn?
- `pp3_min_predictors` — which codes trigger in the first place (see
  [`vignette("acmg-evidence-codes")`](https://ramizkrdnz.github.io/msaVariant/articles/acmg-evidence-codes.md))

Only the last one changes the underlying evidence call. The other two
are presentation.

## Missing data behaves like an off switch

A track is drawn only when its toggle is `TRUE` **and** the bundle
actually carries data for it. A bundle with no REVEL table produces the
same figure as `show_revel = FALSE` — the row is not reserved and then
left blank.

This means you can pass the same explicit toggle set across a whole
cohort of genes with uneven annotation coverage and every figure will
still come out tight.

## Windowing

Not a layer switch, but the other lever for figure density: `window`
sets the residue range on the x-axis, defaulting to ±20 residues around
the variant, clipped to the protein.

``` r

plot_variant_overlay(
  gene          = "DEMO1",
  aligned_fasta = demo_fasta,
  variant_pos   = 21,
  variant_label = "p.R21H",
  window        = c(15, 27),
  show_gnomad   = FALSE,
  show_cadd     = FALSE
)
```

![DEMO1 p.R21H figure restricted to residues
15-27](toggling-layers_files/figure-html/window-1.png)

A narrow window makes individual residue letters legible in a two-column
journal layout; a wide one shows more of the conservation context. Both
stay aligned to the same shared x coordinate system, so every track
lines up with the alignment column above it.

## See also

- [`vignette("annotation-tracks")`](https://ramizkrdnz.github.io/msaVariant/articles/annotation-tracks.md)
  — what each track shows
- [`vignette("colour-schemes")`](https://ramizkrdnz.github.io/msaVariant/articles/colour-schemes.md)
  — restyling what is left after toggling
- [`?plot_variant_overlay`](https://ramizkrdnz.github.io/msaVariant/reference/plot_variant_overlay.md)
  — all arguments
