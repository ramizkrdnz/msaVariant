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
