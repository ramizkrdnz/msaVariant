# msaVariant

`msaVariant` produces publication-quality figures for a **single genetic
variant**: a cross-species protein **multiple sequence alignment** with
per-column **conservation**, stacked **annotation tracks** (ClinVar,
gnomAD, AlphaMissense, REVEL, CADD, protein domains), and an automatic
**ACMG evidence-code strip** (PS1, PM1, PM2, PM5, PP3) computed straight
from the data. It is built on top of
[`ggmsa`](https://github.com/YuLab-SMU/ggmsa) and returns a normal
`ggplot`/`patchwork` object, so every layer, colour, and label stays
fully customisable for figures you can drop into a paper.

## 🔨 Installation

Install the visual base [`ggmsa`](https://github.com/YuLab-SMU/ggmsa)
(Bioconductor), then `msaVariant` from GitHub:

``` r

# ggmsa (required base) + Bioconductor deps
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggmsa")

# msaVariant (development version)
# install.packages("devtools")
devtools::install_github("ramizkrdnz/msaVariant")
```

## 💡 Quick Example

One call assembles the whole figure. The example uses a **local gene
bundle** via
[`import_local_bundle()`](https://ramizkrdnz.github.io/msaVariant/reference/import_local_bundle.md),
so it runs **without any Zenodo download** — point it at the TP53
example bundle shipped in the source repo under
`tests/testthat/fixtures/`, or at any bundle you build with
`data-raw/build_gene_files.R`.

``` r

library(msaVariant)

# 1. Load a local TP53 annotation bundle into the cache (no Zenodo needed)
import_local_bundle("tests/testthat/fixtures/TP53.rds", gene = "TP53")

# 2. Build the figure for the classic pathogenic hotspot TP53 p.R175H
plot_variant_overlay(
  gene          = "TP53",
  aligned_fasta = "tests/testthat/fixtures/tp53_aligned.fasta",
  variant_pos   = 175,
  variant_label = "p.R175H"
)
```

![msaVariant overlay for TP53 p.R175H: cross-species alignment,
conservation, annotation tracks, and ACMG evidence
codes](reference/figures/README-tp53-r175h.png)

For TP53 p.R175H this fires five ACMG codes — **PS1, PM1, PM2, PM5,
PP3** — automatically, from the bundled ClinVar / AlphaMissense / REVEL
/ CADD tables.

## 📚 Learn more

Full articles are on the documentation site,
**<https://ramizkrdnz.github.io/msaVariant/>**:

- [**Get
  Started**](https://ramizkrdnz.github.io/msaVariant/articles/introduction.html)
  — install, load data, and build your first overlay.
- [Full workflow
  tutorial](https://ramizkrdnz.github.io/msaVariant/articles/tutorial_full_workflow.html)
  — end-to-end, from FASTA to finished figure.
- ACMG Evidence Codes — how PS1/PM1/PM2/PM5/PP3 are computed *(in
  preparation)*.
- Colour Schemes — `journal`, `colorblind` (Okabe–Ito), `grayscale`, and
  per-element overrides *(in preparation)*.
- Annotation Tracks — the ClinVar / gnomAD / AlphaMissense / REVEL /
  CADD and domain geoms *(in preparation)*.
- Toggling Layers — show/hide any track; the layout re-flows with no
  gaps *(in preparation)*.
- Bring Your Own Data — overlay your own annotations with
  [`geom_track()`](https://ramizkrdnz.github.io/msaVariant/reference/geom_track.md)
  *(in preparation)*.

## 🧬 Bring-your-own-data mode

No Zenodo bundle required — overlay any per-residue values you already
have:

``` r

my_data <- data.frame(pos = 510:530, score = runif(21))
ggmsa::ggmsa(fa) +
  geom_track(my_data, msa = fa, value = "score", name = "Custom track")
```

## ⚠️ Data upload pending

The Zenodo data deposit is **not yet published**, so remote fetches
([`fetch_gene_data()`](https://ramizkrdnz.github.io/msaVariant/reference/fetch_gene_data.md)
and the gene-symbol annotation layers) will currently fail. Until the
deposit is live, use your own data via *Bring-your-own-data mode* above,
or load a local bundle with
[`import_local_bundle()`](https://ramizkrdnz.github.io/msaVariant/reference/import_local_bundle.md)
(as in the Quick Example). This note will be removed once the deposit is
uploaded and wired in.

## 📄 License

`msaVariant` is released under the **Artistic-2.0** license.

Note that data fetched via
[`geom_alphamissense()`](https://ramizkrdnz.github.io/msaVariant/reference/geom_alphamissense.md)
is licensed **CC-BY-NC-SA 4.0**, which restricts commercial use of works
derived from it. Check the license of every annotation source you
display before reusing a figure.

## 🏃 Author

*Author and maintainer details are pending and will be added here.*
