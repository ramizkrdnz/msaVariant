# Resolve the full color specification for a variant-overlay figure

Resolve the full color specification for a variant-overlay figure

## Usage

``` r
resolve_plot_colors(
  color_scheme = "journal",
  acmg_colors = NULL,
  variant_highlight_color = NULL,
  track_palettes = NULL,
  msa_color_scheme = NULL
)
```

## Arguments

- color_scheme:

  One of \`"journal"\` (default), \`"colorblind"\`, \`"grayscale"\`.

- acmg_colors:

  Optional named character vector overriding ACMG chip colors (names
  among PS1/PM1/PM2/PM5/PP3). Partial is allowed.

- variant_highlight_color:

  Optional single color for the variant highlight fill.

- track_palettes:

  Optional named list with any of \`clinvar\`, \`gnomad\`, \`cadd\`,
  \`alphamissense\`, \`revel\`.

- msa_color_scheme:

  Optional MSA palette: a scheme name (\`"journal"\`, \`"colorblind"\`,
  \`"clustal"\`, \`"grayscale"\`) or a named per-residue color vector.

## Value

A named list of resolved colors.
