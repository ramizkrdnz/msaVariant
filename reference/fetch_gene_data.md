# Fetch the combined annotation file for a gene

Downloads (or reads from local cache) the per-gene annotation bundle
from the Zenodo deposit. Returns the deserialized list as described in
\`DATA_FORMAT_SPEC.md\`.

## Usage

``` r
fetch_gene_data(
  gene,
  force_refresh = FALSE,
  validate = TRUE,
  verify_checksum = TRUE,
  quiet = FALSE
)
```

## Arguments

- gene:

  HGNC gene symbol (e.g. \`"PATL1"\`).

- force_refresh:

  If \`TRUE\`, redownload even if cached.

- validate:

  If \`TRUE\` (default), validate the file against the package's format
  spec before returning.

- verify_checksum:

  If \`TRUE\` (default), and a local \`MANIFEST.tsv\` is present, verify
  the bundle's sha256 against the manifest entry; a mismatched file is
  treated as corrupt, removed from the cache, and \`NULL\` is returned.
  Has no effect when no manifest is present.

- quiet:

  If \`TRUE\`, suppress "Downloading..." messages.

## Value

A named list with 7 elements (\`meta\`, \`domains\`, \`clinvar\`,
\`gnomad\`, \`alphamissense\`, \`revel\`, \`cadd\`), or \`NULL\` with a
warning on failure.

## Details

Most users will not call this directly; the \`get_domains()\`,
\`get_clinvar()\`, etc. helpers and the \`geom\_\*()\` layers route
through it transparently.

## Examples

``` r
## Runnable with the shipped synthetic DEMO1 bundle (no network).
## A temporary cache keeps the example off your real cache directory.
Sys.setenv(MSAVARIANT_CACHE = tempfile("msaVariant_cache_"))
import_local_bundle(
    system.file("extdata", "DEMO1.rds", package = "msaVariant"),
    gene = "DEMO1"
)
#> Cache directory does not exist yet (no annotations have been downloaded).
#> Imported DEMO1 bundle -> /tmp/Rtmp9JzdHU/msaVariant_cache_1a73584ba281/0.1.0/DEMO1.rds
b <- fetch_gene_data("DEMO1")
names(b)
#> [1] "meta"          "domains"       "clinvar"       "gnomad"       
#> [5] "alphamissense" "revel"         "cadd"         

if (FALSE) { # \dontrun{
## Real genes are downloaded from the Zenodo data deposit.
b <- fetch_gene_data("TP53")
} # }
```
