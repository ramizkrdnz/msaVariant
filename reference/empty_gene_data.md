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
