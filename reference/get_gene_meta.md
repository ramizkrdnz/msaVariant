# Get gene metadata

Returns the \`meta\` slice — one row with UniProt ID, protein length,
Ensembl IDs, build date, and source-version stamps.

## Usage

``` r
get_gene_meta(gene, force_refresh = FALSE)
```

## Arguments

- gene:

  HGNC gene symbol.

- force_refresh:

  Redownload even if cached.

## Value

A 1-row \`data.frame\`; returns \`NULL\` on failure.

## Examples

``` r
Sys.setenv(MSAVARIANT_CACHE = tempfile("msaVariant_cache_"))
import_local_bundle(
    system.file("extdata", "DEMO1.rds", package = "msaVariant"),
    gene = "DEMO1"
)
#> Cache directory does not exist yet (no annotations have been downloaded).
#> Imported DEMO1 bundle -> /tmp/Rtmp9JzdHU/msaVariant_cache_1a73339ec666/0.1.0/DEMO1.rds
get_gene_meta("DEMO1")
#>    gene uniprot_id protein_length     ensembl_gene_id ensembl_transcript_id
#> 1 DEMO1  DEMO00001             40 ENSGDEMO00000000001   ENSTDEMO00000000001
#>   refseq_protein   mane_select build_date
#> 1  NP_DEMO0001.1 NM_DEMO0001.1 2026-07-08
#>                                                             source_versions
#> 1 SYNTHETIC — illustrative only, not real AlphaMissense/CADD/ClinVar/gnomAD
```
