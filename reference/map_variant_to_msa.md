# Map a vector of variant residue positions to MSA columns

Thin wrapper around \`build_msa_coord_map()\` that vectorises the lookup
and warns when a position is out of range.

## Usage

``` r
map_variant_to_msa(positions, msa, ref_name = NULL)
```

## Arguments

- positions:

  Integer vector of reference protein positions.

- msa, ref_name:

  See \`build_msa_coord_map()\`.

## Value

Integer vector of MSA column indices (NA where mapping fails).

## Examples

``` r
fa <- system.file("extdata", "demo_aligned.fasta", package = "msaVariant")
map_variant_to_msa(c(21, 30), fa, ref_name = "DEMO1_HUMAN")
#> [1] 21 30
```
