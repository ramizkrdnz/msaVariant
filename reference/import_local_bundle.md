# Import a locally-built gene bundle into the cache

Copies the supplied \`.rds\` file into the package's cache directory
(see \[cache_location()\]) so that \[fetch_gene_data()\] and
\[plot_variant_overlay()\] can find it.

## Usage

``` r
import_local_bundle(
  path,
  gene = NULL,
  overwrite = TRUE,
  verify_checksum = TRUE
)
```

## Arguments

- path:

  Path to a \`.rds\` file produced by your build script. The file must
  conform to the gene-data schema enforced by \[validate_gene_data()\].

- gene:

  Optional gene symbol; if \`NULL\` (default) the symbol is read from
  \`bundle\$meta\$gene\` inside the file. Used as the destination
  filename: \`\<GENE\>.rds\`.

- overwrite:

  Logical, default \`TRUE\`. If \`FALSE\` and a bundle for this gene is
  already present, an error is thrown.

- verify_checksum:

  Logical, default \`TRUE\`. If a local \`MANIFEST.tsv\` is present and
  lists this gene, the imported file's sha256 is compared against the
  manifest entry and a warning is emitted on mismatch. The file is still
  imported (you asked for it explicitly); the warning simply flags a
  possible corruption or a stale manifest. No effect when no manifest is
  present.

## Value

Invisibly, the destination path of the imported bundle.

## Examples

``` r
## Import the shipped synthetic DEMO1 bundle into a temporary cache.
Sys.setenv(MSAVARIANT_CACHE = tempfile("msaVariant_cache_"))
import_local_bundle(
  system.file("extdata", "DEMO1.rds", package = "msaVariant"),
  gene = "DEMO1"
)
#> Cache directory does not exist yet (no annotations have been downloaded).
#> Imported DEMO1 bundle -> /tmp/RtmpUwZ6Qa/msaVariant_cache_1a0c3d39973c/0.1.0/DEMO1.rds
fetch_gene_data("DEMO1")$meta
#>    gene uniprot_id protein_length     ensembl_gene_id ensembl_transcript_id
#> 1 DEMO1  DEMO00001             40 ENSGDEMO00000000001   ENSTDEMO00000000001
#>   refseq_protein   mane_select build_date
#> 1  NP_DEMO0001.1 NM_DEMO0001.1 2026-07-08
#>                                                             source_versions
#> 1 SYNTHETIC — illustrative only, not real AlphaMissense/CADD/ClinVar/gnomAD
```
