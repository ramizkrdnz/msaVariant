# Default colour scale for ACMG pathogenicity tiers

Provides a colourblind-safe 5-tier palette ordered from pathogenic (red)
to benign (blue), with grey for VUS.

## Usage

``` r
scale_fill_pathogenicity(...)
```

## Arguments

- ...:

  Passed to \`ggplot2::scale_fill_manual()\`.

## Value

A \`ScaleDiscrete\` object.
