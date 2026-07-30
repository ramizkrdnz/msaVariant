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

## Examples

``` r
library(ggplot2)
d <- data.frame(x = c("Pathogenic", "VUS", "Benign"), y = 1)
ggplot(d, aes(x, y, fill = x)) +
    geom_col() +
    scale_fill_pathogenicity()
```
