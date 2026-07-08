# Tutorial: From raw protein sequences to an annotated MSA figure

## What this tutorial covers

The complete workflow:

1.  Obtain protein orthologs from UniProt (unaligned).
2.  Run a multiple sequence alignment in R.
3.  Visualise with `msaVariant`, including the patient’s variant and
    database-derived annotation tracks.

We use the **PATL1** human gene and seven orthologs as the running
example, because that’s the case `msaVariant` was originally written
for. The workflow generalises to any single-protein investigation.

## Step 1 — Install the alignment software

`msaVariant` does not perform alignment itself. Choose any mainstream
protein-MSA tool; we recommend `msa` (Bioconductor, ships
MUSCLE/ClustalW/ClustalOmega) because it stays inside R:

``` r

if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("msa")
```

Standalone alternatives (run from the shell, then load with
[`Biostrings::readAAMultipleAlignment`](https://rdrr.io/pkg/Biostrings/man/MultipleAlignment-class.html)):

- `mafft` — `apt-get install mafft` / `brew install mafft`
- `muscle` — `apt-get install muscle`
- `clustalo` — `apt-get install clustalo`

For this tutorial we use `msa::msa()` with MUSCLE.

## Step 2 — Get the unaligned ortholog sequences

You need a FASTA file with the protein sequences of your gene of
interest in human plus several orthologs spanning the evolutionary
distance you want to argue about. The faster and more reliable route is
to pull canonical sequences from UniProt by accession:

``` r

library(Biostrings)

# UniProt accessions for PATL1 across eight species
accessions <- c(
  "PATL1_HUMAN" = "Q86TB9",        # Homo sapiens
  "PATL1_PANTR" = "A0A2I3RIF4",    # Pan troglodytes
  "PATL1_MOUSE" = "Q3TC46",        # Mus musculus
  "PATL1_RAT"   = "B5DF93",        # Rattus norvegicus
  "PATL1_XENLA" = "Q32N92",        # Xenopus laevis
  "PATL1_DANRE" = "A2RRV3",        # Danio rerio
  "PATR1_DROME" = "Q9VEN9",        # Drosophila melanogaster
  "PATR1_CAEEL" = "Q20374"         # Caenorhabditis elegans
)

# Pull each FASTA from UniProt. Each download is ~1 KB.
seqs <- lapply(accessions, function(acc) {
  url <- sprintf("https://rest.uniprot.org/uniprotkb/%s.fasta", acc)
  s   <- readAAStringSet(url)
  s[1]   # first (canonical) sequence in the file
})

# Rename to readable IDs and concatenate
unaligned <- do.call(c, seqs)
names(unaligned) <- names(accessions)

# Save the unaligned FASTA in case you want to inspect it
writeXStringSet(unaligned, "patl1_unaligned.fasta")
```

> **Tip**: If you already have a FASTA file from your collaborators,
> skip the download and just `readAAStringSet("yourfile.fasta")`.

## Step 3 — Align them

``` r

library(msa)

aln <- msa(unaligned, method = "Muscle", type = "protein")
aln

# Convert to a plain FASTA file that `msaVariant` can read
aa_aligned <- as(aln, "AAMultipleAlignment")
writeXStringSet(as(aa_aligned, "AAStringSet"),
                filepath = "patl1_aligned.fasta")
```

You should see a multiple alignment of length 800-1300 columns (varies
with the alignment software’s gap insertions).

## Step 4 — Pick a residue window of interest

The full alignment is too wide to render as a single figure. Pick a
residue range in the human reference where your variant sits. For the
PATL1 K518fs case, we focus on residues 490..550:

``` r

library(msaVariant)

# Use the package's coord-map utility to convert
# "human residue 490..550" into "alignment columns N..M"
map <- build_msa_coord_map("patl1_aligned.fasta",
                            ref_name = "PATL1_HUMAN")
start_col <- map$msa_col[match(490, map$residue_pos)]
end_col   <- map$msa_col[match(550, map$residue_pos)]
c(start_col, end_col)
```

Extract that window into a smaller FASTA for plotting:

``` r

aa <- readAAStringSet("patl1_aligned.fasta")
window <- subseq(aa, start = start_col, end = end_col)
writeXStringSet(window, "patl1_window.fasta")
```

## Step 5 — Define the patient’s variant

``` r

patient <- data.frame(
  pos         = 518 - 490 + 1,    # position within the window
  pos_end     = 550 - 490 + 1,    # frameshift extends to end of window
  label       = "K518fs",
  consequence = factor("frameshift",
                       levels = c("frameshift","missense","nonsense")),
  stringsAsFactors = FALSE
)
```

The `pos` column is in **window-local** residue numbering (1 = first
residue in the window). The package translates that to alignment columns
automatically.

## Step 6 — The figure

``` r

library(ggmsa)

ggmsa("patl1_window.fasta",
      char_width = 0.6, seq_name = TRUE) +
  geom_variant(patient,                    msa = "patl1_window.fasta") +
  geom_domain(gene = "PATL1",              msa = "patl1_window.fasta") +
  geom_clinvar(gene = "PATL1",             msa = "patl1_window.fasta") +
  geom_alphamissense(gene = "PATL1",       msa = "patl1_window.fasta") +
  geom_gnomad(gene = "PATL1",              msa = "patl1_window.fasta")
```

The first time you call each `geom_*(gene = ...)`, the package downloads
the per-gene annotation slice from the Zenodo deposit (~50–200 KB each)
and caches it. Subsequent runs are instant.

## Common pitfalls

**“Residue position X not found”** — your variant is outside the window
you extracted. Either widen the window or check your residue numbering
(full-protein vs window-local).

**Sequence names don’t match** — `geom_*` defaults to the first sequence
in the FASTA as the reference. If your human PATL1 isn’t first, supply
`ref_name = "PATL1_HUMAN"` to every layer.

**Empty annotation track** — `geom_gnomad(gene = "FOO")` returned
`NULL`. Either the gene isn’t in our Zenodo deposit, or your
internet/firewall blocked the download. Try
[`clear_cache(); get_gnomad("FOO")`](https://ramizkrdnz.github.io/msaVariant/reference/clear_cache.md)
to see the actual error.

**The figure is dominated by the gnomAD/CADD/AlphaMissense tracks** —
they each take ~1 row of vertical space. For 8 sequences + 5 annotation
tracks, set `fig.height = 7` or larger.

## Adjusting for non-PATL1 use cases

For any other gene, the only changes are:

1.  The UniProt accessions in Step 2.
2.  The gene symbol in `geom_*(gene = "...")` calls.
3.  The residue window in Steps 4-5.

Everything else stays the same.
