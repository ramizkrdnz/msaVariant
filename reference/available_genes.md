# List the genes for which annotation data is available

Returns the set of gene symbols the package can load. In the grouped
deposit layout (the default), the authoritative source is the
gene-\>group index shipped in the package
(\`inst/extdata/gene_group_index.tsv\`): it lists every gene and needs
no download or prior caching. In the legacy one-file-per-gene layout
(grouped mode disabled, or no index present), the local \`MANIFEST.tsv\`
in the cache is used instead, listing the genes whose bundles have been
downloaded or imported.

## Usage

``` r
available_genes()
```

## Value

A sorted character vector of HGNC gene symbols. Empty if no index and no
manifest are found.

## Details

When neither source is available, an empty character vector is returned.

## See also

\[fetch_gene_data()\], \[import_local_bundle()\], \[cache_summary()\]

## Examples

``` r
## Grouped mode (default): reads the gene->group index shipped in the
## package, so this works offline, with no cached data.
available_genes()
#> [1] "DEMO1"  "SYNTH1" "SYNTH2" "TP53"   "WDR31" 
```
