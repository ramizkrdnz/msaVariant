# List genes available in the local data manifest

Reads the \`MANIFEST.tsv\` that ships alongside the cached per-gene
bundles (see \[cache_location()\]) and returns the gene symbols it
lists. This is the set of genes for which annotation data can be loaded
without a further download.

## Usage

``` r
available_genes()
```

## Value

A sorted character vector of HGNC gene symbols. Empty if no manifest is
found.

## Details

If no manifest is present (for example before any data has been
downloaded or imported), an empty character vector is returned.

## See also

\[fetch_gene_data()\], \[import_local_bundle()\], \[cache_summary()\]

## Examples

``` r
if (FALSE) { # \dontrun{
available_genes()
} # }
```
