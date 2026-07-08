# msaVariant

> Clinical-genetics multiple sequence alignment visualisation with
> variant overlay. v0.2.0 — Kaplan Lab.

`msaVariant` extends [`ggmsa`](https://github.com/YuLab-SMU/ggmsa)
with annotation layers tailored to clinical and rare-disease
genetics: variants, protein domains, gnomAD allele frequencies,
ClinVar significance, AlphaMissense/REVEL/CADD pathogenicity
scores.

## Installation

```r
# Install ggmsa first (required for the visual base)
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("ggmsa")

# Then msaVariant
# (Bioconductor route once published:)
BiocManager::install("msaVariant")
# (GitHub route in the meantime:)
# devtools::install_github("KaplanLab/msaVariant")
```

## Quick start

```r
library(ggmsa)
library(msaVariant)

fa <- "patl1_orthologs.fasta"  # your MSA

ggmsa(fa, char_width = 0.6, seq_name = TRUE) +
  geom_variant(data.frame(pos = 518, label = "K518fs"),
               msa = fa) +
  geom_domain(gene = "PATL1",         msa = fa) +
  geom_clinvar(gene = "PATL1",        msa = fa) +
  geom_alphamissense(gene = "PATL1",  msa = fa) +
  geom_gnomad(gene = "PATL1",         msa = fa)
```

Five layers, one gene symbol per annotation, no manual data
download. The package fetches per-gene data from a Zenodo deposit
the first time you use a gene and caches it locally.

> **⚠️ Data upload pending.** The Zenodo data deposit is not yet
> published, so remote fetches (`fetch_gene_data()` and the
> gene-symbol annotation layers above) will currently fail. Until
> the deposit is live you can still use every feature with your own
> data via *Bring-your-own-data mode* (below) or by importing a
> local bundle with `import_local_bundle()`. This note will be
> removed once the deposit is uploaded and wired in.

## Data architecture

| Annotation       | Where it comes from                              |
|------------------|--------------------------------------------------|
| Conservation     | Computed from your MSA                           |
| Domains          | Zenodo deposit (InterPro/Pfam slice)             |
| ClinVar          | Zenodo deposit (NCBI snapshot)                   |
| gnomAD           | Zenodo deposit (per-residue summary, v4.1)       |
| AlphaMissense    | Zenodo deposit (per-residue mean/max scores)     |
| REVEL            | Zenodo deposit                                   |
| CADD             | Zenodo deposit                                   |

The data deposit is versioned by the Kaplan Lab on Zenodo. When
gnomAD releases a new version, we re-run our build scripts
(`data-raw/build_zenodo_deposit.R`), upload a new deposit version,
and release a new `msaVariant` package version that points to it.

## Bring-your-own-data mode

If you have proprietary or unpublished annotation:

```r
my_data <- data.frame(pos = 510:530, score = runif(21))
ggmsa(fa) +
  geom_track(my_data, msa = fa, value = "score",
             name = "Custom track")
```

## Citation

If you use this package in published work please cite:

* `msaVariant` (this package) — Kaplan Lab.
* `ggmsa` — Zhou L et al. Briefings in Bioinformatics 23(4):bbac222 (2022).
* The data sources you used — see the Zenodo deposit README for
  individual citations.

## License

Artistic-2.0. Note that data fetched via `geom_alphamissense()` is
licensed CC-BY-NC-SA 4.0, which restricts commercial use of works
derived from it.

## For maintainers

To rebuild the data deposit (gnomAD release, ClinVar refresh,
etc.), read `data-raw/ZENODO_UPLOAD.md`.
