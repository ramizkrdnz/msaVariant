# Construct an empty per-gene object conforming to the spec

Builds a placeholder list with all 7 required tables present but empty
(zero rows, correct columns and types). Useful as a template when
writing custom annotation files or as a fallback.

## Usage

``` r
empty_gene_data(
  gene = "UNKNOWN",
  uniprot_id = NA_character_,
  protein_length = 0L
)
```

## Arguments

- gene:

  HGNC symbol for the \`meta\$gene\` field.

- uniprot_id:

  UniProt accession.

- protein_length:

  Integer length.

## Value

A list conforming to \`DATA_FORMAT_SPEC.md\`.

## Examples

``` r
skeleton <- empty_gene_data("MYGENE",
    uniprot_id = "P00000",
    protein_length = 393L
)
lapply(skeleton, names)
#> $meta
#> [1] "gene"                  "uniprot_id"            "protein_length"       
#> [4] "ensembl_gene_id"       "ensembl_transcript_id" "build_date"           
#> [7] "source_versions"      
#> 
#> $domains
#> [1] "start"     "end"       "name"      "accession" "source"   
#> 
#> $clinvar
#> [1] "pos"           "aa_ref"        "aa_alt"        "aa_change"    
#> [5] "significance"  "review_status" "clinvar_id"   
#> 
#> $gnomad
#> [1] "pos"         "aa_ref"      "aa_alt"      "aa_change"   "consequence"
#> [6] "af_joint"    "ac_joint"    "an_joint"    "filter"     
#> 
#> $alphamissense
#> [1] "pos"       "aa_ref"    "aa_alt"    "aa_change" "am_score"  "am_class" 
#> 
#> $revel
#> [1] "pos"         "aa_ref"      "aa_alt"      "aa_change"   "revel_score"
#> 
#> $cadd
#> [1] "pos"         "aa_ref"      "aa_alt"      "aa_change"   "consequence"
#> [6] "cadd_raw"    "cadd_phred" 
#> 
```
