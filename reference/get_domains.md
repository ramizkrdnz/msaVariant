# Get InterPro/Pfam domains for a gene

Returns the \`domains\` slice of the per-gene Zenodo bundle.

## Usage

``` r
get_domains(gene, force_refresh = FALSE)
```

## Arguments

- gene:

  HGNC gene symbol (e.g. \`"PATL1"\`).

- force_refresh:

  Redownload even if cached.

## Value

A \`data.frame\` with columns \`start\`, \`end\`, \`name\`,
\`accession\`, \`source\` (factor), and optionally \`evidence\`. Returns
\`NULL\` if download failed.

## Examples

``` r
## Runnable with the shipped synthetic DEMO1 bundle (no network).
Sys.setenv(MSAVARIANT_CACHE = tempfile("msaVariant_cache_"))
import_local_bundle(
    system.file("extdata", "DEMO1.rds", package = "msaVariant"),
    gene = "DEMO1"
)
#> Cache directory does not exist yet (no annotations have been downloaded).
#> Imported DEMO1 bundle -> /tmp/Rtmp9JzdHU/msaVariant_cache_1a733468ecda/0.1.0/DEMO1.rds
get_domains("DEMO1")
#>   start end                   name  accession   source
#> 1     5  35 Demo functional domain PFDEMO0001     Pfam
#> 2     2  10 Demo N-terminal region IPRDEMO001 InterPro
```
