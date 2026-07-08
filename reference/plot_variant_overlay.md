# Generate a publication-grade variant-overlay figure

Composes a multi-panel figure (variant header, MSA, annotation tracks,
legend dashboard) for a clinical variant in the context of evolutionary
conservation and population/clinical annotations. The gene's annotation
data is read from the package cache via \[fetch_gene_data()\]. The user
supplies only the gene symbol, the aligned FASTA, and the variant
position.

## Usage

``` r
plot_variant_overlay(
  gene,
  aligned_fasta,
  variant_pos,
  variant_label = NULL,
  window = NULL,
  ref_name = NULL,
  title = NULL,
  subtitle = NULL,
  aa_change = NULL,
  pp3_min_predictors = 2L,
  color_scheme = "journal",
  acmg_colors = NULL,
  variant_highlight_color = NULL,
  track_palettes = NULL,
  msa_color_scheme = NULL,
  show_acmg = TRUE,
  show_domain = TRUE,
  show_clinvar = TRUE,
  show_gnomad = TRUE,
  show_cadd = TRUE,
  show_alphamissense = TRUE,
  show_revel = TRUE,
  acmg_codes = c("PS1", "PM1", "PM2", "PM5", "PP3")
)
```

## Arguments

- gene:

  Gene symbol (string). Must match the name of a \`.rds\` bundle in the
  cache. Use \[import_local_bundle()\] to add a local bundle, or
  \[fetch_gene_data()\] to download from Zenodo.

- aligned_fasta:

  Path to an aligned FASTA file (e.g. MUSCLE output). The reference
  sequence (matching the canonical UniProt isoform) should be one of the
  sequences in the file.

- variant_pos:

  Residue position in the reference (UniProt canonical numbering). Can
  also be a string like \`"R175H"\` or \`"p.R175H"\` from which the
  position is parsed.

- variant_label:

  Optional display label (e.g. \`"p.R175H"\`). If \`NULL\`, derived from
  \`variant_pos\` and the reference residue.

- window:

  Length-2 integer vector \`c(start, end)\` defining the residue range
  to display. Defaults to ±20 residues around \`variant_pos\`, clipped
  to the protein.

- ref_name:

  Name of the reference sequence in the FASTA. If \`NULL\`, the sequence
  whose name contains the gene symbol is used; otherwise the first
  sequence.

- title:

  Plot title.

- subtitle:

  Plot subtitle.

- aa_change:

  Optional variant identifier used for the ACMG evidence strip, e.g.
  \`"R175H"\`. If \`NULL\`, it is parsed from \`variant_label\` when
  that carries a full ref+pos+alt substitution. ACMG codes that require
  the alternate allele are only shown when a complete \`aa_change\` is
  available.

- pp3_min_predictors:

  Integer 1–3 (default 2) forwarded to \[compute_acmg_codes()\]: PP3
  fires when at least this many of the three PP3 predictors
  (AlphaMissense/REVEL/CADD) pass their threshold.

- color_scheme:

  Figure color preset: \`"journal"\` (default, publication-quality,
  colorblind-safe, grayscale-legible, muted), \`"colorblind"\` (explicit
  Okabe-Ito), or \`"grayscale"\` (pure black-and-white for print).

- acmg_colors:

  Optional named vector overriding ACMG chip colors (names among
  \`PS1\`,\`PM1\`,\`PM2\`,\`PM5\`,\`PP3\`); partial allowed. \`NULL\`
  uses the scheme.

- variant_highlight_color:

  Optional single color for the variant highlight fill; \`NULL\` uses
  the scheme.

- track_palettes:

  Optional named list overriding individual track palettes: any of
  \`clinvar\`, \`gnomad\`, \`cadd\`, \`alphamissense\`, \`revel\`.
  \`NULL\` uses the scheme.

- msa_color_scheme:

  Optional MSA coloring: a scheme name (\`"journal"\`, \`"colorblind"\`,
  \`"clustal"\`, \`"grayscale"\`) or an explicit named per-residue color
  vector. \`NULL\` uses the scheme.

- show_acmg, show_domain, show_clinvar, show_gnomad, show_cadd,
  show_alphamissense, show_revel:

  Logical toggles (all \`TRUE\` by default) for each annotation layer.
  Setting one to \`FALSE\` fully removes that layer and the figure
  re-flows (patchwork heights and row assignments adjust to the number
  of active layers — no blank gaps). The variant header and the MSA are
  always shown.

- acmg_codes:

  Character vector of ACMG codes to \*display\* (default all five:
  \`c("PS1","PM1","PM2","PM5","PP3")\`). The full ACMG computation
  always runs internally; only codes in this set that were actually
  triggered are drawn as chips.

## Value

A \`patchwork\` plot object. Save with \`ggplot2::ggsave()\`.

## Examples

``` r
## Runnable with the synthetic DEMO1 example bundle shipped in the
## package (fabricated data; not real predictions). Use a temporary
## cache so the example does not touch your real cache directory.
Sys.setenv(MSAVARIANT_CACHE = tempfile("msaVariant_cache_"))
import_local_bundle(
  system.file("extdata", "DEMO1.rds", package = "msaVariant"),
  gene = "DEMO1"
)
#> Cache directory does not exist yet (no annotations have been downloaded).
#> Imported DEMO1 bundle -> /tmp/RtmpB0BDjI/msaVariant_cache_1a583882f1df/0.1.0/DEMO1.rds
p <- plot_variant_overlay(
  gene          = "DEMO1",
  aligned_fasta = system.file("extdata", "demo_aligned.fasta",
                              package = "msaVariant"),
  variant_pos   = 21,
  variant_label = "p.R21H"
)

if (FALSE) { # \dontrun{
## Save a publication-resolution figure
ggplot2::ggsave("demo1_R21H.png", p,
                width = 15, height = 9, dpi = 300, bg = "white")
} # }
```
