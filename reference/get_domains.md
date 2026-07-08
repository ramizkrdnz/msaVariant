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
