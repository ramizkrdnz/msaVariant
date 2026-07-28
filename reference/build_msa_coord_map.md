# Build a residue -\> MSA column lookup

For one row of a multiple sequence alignment (the "reference" sequence),
build a map from ungapped residue position (1, 2, 3, ...) to alignment
column.

## Usage

``` r
build_msa_coord_map(msa, ref_name = NULL)
```

## Arguments

- msa:

  An object understood by ggmsa: an AAMultipleAlignment,
  DNAMultipleAlignment, AAStringSet, DNAStringSet, character vector of
  equal-length sequences, or a path to a FASTA file.

- ref_name:

  Name of the reference sequence in the MSA. The default \`NULL\` uses
  the first sequence.

## Value

A data.frame with columns \`residue_pos\`, \`aa\`, \`msa_col\`. Gap
columns in the reference are skipped.

## Examples

``` r
fa <- system.file("extdata", "demo_aligned.fasta", package = "msaVariant")
cmap <- build_msa_coord_map(fa, ref_name = "DEMO1_HUMAN")
head(cmap)
#>   residue_pos aa msa_col
#> 1           1  M       1
#> 2           2  K       2
#> 3           3  T       3
#> 4           4  A       4
#> 5           5  Y       5
#> 6           6  I       6
```
