# Introduction to msaVariant

## What this package does

`msaVariant` overlays clinical-genetics evidence onto a multiple
sequence alignment (MSA) in a single ggplot2 figure. Each layer paints a
different kind of evidence:

- **[`geom_variant()`](https://ramizkrdnz.github.io/msaVariant/reference/geom_variant.md)**
  — the patient’s variant(s)
- **[`geom_domain()`](https://ramizkrdnz.github.io/msaVariant/reference/geom_domain.md)**
  — InterPro/Pfam protein domains
- **[`geom_gnomad()`](https://ramizkrdnz.github.io/msaVariant/reference/geom_gnomad.md)**
  — gnomAD allele frequencies
- **[`geom_clinvar()`](https://ramizkrdnz.github.io/msaVariant/reference/geom_clinvar.md)**
  — ClinVar pathogenicity calls
- **[`geom_alphamissense()`](https://ramizkrdnz.github.io/msaVariant/reference/geom_alphamissense.md)**
  — AlphaMissense in silico scores
- **[`geom_revel()`](https://ramizkrdnz.github.io/msaVariant/reference/geom_revel.md)
  /
  [`geom_cadd()`](https://ramizkrdnz.github.io/msaVariant/reference/geom_cadd.md)**
  — REVEL / CADD scores
- **[`geom_track()`](https://ramizkrdnz.github.io/msaVariant/reference/geom_track.md)**
  — generic per-residue track (bring your own data)

## How the data layer works

The five annotation-database layers (`geom_gnomad`, `geom_clinvar`,
`geom_alphamissense`, `geom_revel`, `geom_cadd`) and `geom_domain` fetch
their data from a Zenodo deposit the first time you use a gene, and
cache it locally for next time. The package code itself ships with no
clinical-genetics data — only the visualisation machinery.

This means:

1.  The package install is small and fast.
2.  You always get the data version we pinned to the package release, so
    the figure is reproducible.
3.  Updating the data (when gnomAD releases v5, for example) is a matter
    of running our build scripts and re-uploading to Zenodo.

The first time you call e.g. `geom_gnomad(gene = "PATL1", msa = fa)`,
the package will download ~100 KB from Zenodo and cache it at
`tools::R_user_dir("msaVariant", "cache")`. Subsequent calls are
instant.

## The PATL1 worked example

``` r

library(ggmsa)
library(msaVariant)

# 1. Your alignment: PATL1 orthologs around residue 518
fa <- system.file("extdata", "patl1_orthologs.fasta",
                  package = "msaVariant")

# 2. The patient's variant
patient <- data.frame(
  pos         = 29,                # K518 in the demo stretch
  pos_end     = 61,                # frameshift through end of window
  label       = "K518fs",
  consequence = factor("frameshift",
                       levels = c("frameshift","missense","nonsense"))
)

# 3. The figure
ggmsa(fa, char_width = 0.6, seq_name = TRUE) +
  geom_variant(patient,                  msa = fa) +
  geom_domain(gene = "PATL1",            msa = fa) +
  geom_clinvar(gene = "PATL1",           msa = fa) +
  geom_alphamissense(gene = "PATL1",     msa = fa) +
  geom_gnomad(gene = "PATL1",            msa = fa)
```

That’s it. Five layers, one gene symbol per annotation, no manual data
download.

## Bringing your own data

If you already have a custom annotation (a lab-internal database, a
paper’s supplementary table, etc.), pass it via `data =` instead of
`gene =`:

``` r

my_annotations <- read.csv("my_annotations.csv")
ggmsa(fa) +
  geom_track(my_annotations,
             msa = fa, value = "score",
             name = "My score")
```

The `data =` and `gene =` modes are mutually exclusive — pass one or the
other, never both.

## Cache management

``` r

cache_location()       # where the cache lives
clear_cache()          # nuke everything
clear_cache("gnomad")  # only the gnomAD slice
```

## Citation

If you use `msaVariant` in published work, please cite both the package
and the underlying data sources (see the Zenodo deposit README for the
full citation list).
