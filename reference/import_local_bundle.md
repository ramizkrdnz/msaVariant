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
if (FALSE) { # \dontrun{
import_local_bundle("/path/to/TP53.rds")
import_local_bundle("/path/to/WDR31.rds")
p <- plot_variant_overlay(
  gene = "TP53",
  aligned_fasta = "tp53_aligned.fasta",
  variant_pos = 175
)
} # }
```
