# Validate a per-gene annotation object against the package spec

Checks that the object conforms to the structure defined in
\`DATA_FORMAT_SPEC.md\`: presence of all 7 named elements, required
columns, correct types, and no \`NA\` in required columns. Returns a
list with \`valid\` (logical) and \`message\` (character) describing all
issues found.

## Usage

``` r
validate_gene_data(obj, strict = FALSE)
```

## Arguments

- obj:

  The deserialized object from a per-gene \`.rds\` file.

- strict:

  If \`TRUE\`, throw an error on validation failure instead of returning
  a list. Default \`FALSE\`.

## Value

If \`strict = FALSE\`: a list with elements \`valid\` (logical) and
\`message\` (character vector of issues, or NULL if valid). If \`strict
= TRUE\`: invisibly \`TRUE\` on success, or \`stop()\` with a formatted
message on failure.

## Details

Used internally by \`fetch_gene_data()\` after every cache hit or
download. Also useful for the lab's build pipeline to validate files
before uploading to Zenodo.
