# msaVariant — Zenodo deposit instructions

This walks through building and publishing the annotation-data
deposit that the `msaVariant` package consumes.

Most days, you don't do this — users `install.packages("msaVariant")`
and per-gene files are downloaded on demand. You only run this:

- Once, for the initial v0.1 release.
- Quarterly (or when one of the underlying databases issues a
  significant update worth refreshing for).

---

## Prerequisites

| Need | Detail |
|---|---|
| Hardware | Linux/macOS, ~500 GB free disk, ≥ 16 GB RAM |
| Time | 6–12 hours mostly unattended |
| Software | R ≥ 4.1, plus `BiocManager::install("Biostrings")` and `install.packages(c("readr","httr","jsonlite","digest","dplyr","stringr"))` |
| Network | Stable broadband; the CADD download alone is ~300 GB |
| Zenodo account | Free at https://zenodo.org/signup/ (use ORCID-linked institutional account if available) |

---

## Step 1 — Verify source URLs before starting

Open `data-raw/build_scripts/build_*.R` and check the URLs at the
top of each script still resolve. URLs that may shift between
releases are flagged with `URL:` comments. If a source has
reorganized, update the URL before running the master script.

Quick health-check loop:

```bash
for url in \
  https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/protein2ipr.dat.gz \
  https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz \
  https://zenodo.org/record/8208688/files/AlphaMissense_aa_substitutions.tsv.gz \
  https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip \
  https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs_inclAnno.tsv.gz \
  https://rest.uniprot.org/uniprotkb/stream?query=organism_id:9606+AND+reviewed:true
do
  echo "Checking $url"
  curl -sI -o /dev/null -w "%{http_code}\n" -L "$url"
done
```

Any 404s — fix in the corresponding `build_*.R` script before
proceeding.

---

## Step 2 — Run the build

From the package root:

```bash
cd path/to/msaVariant
Rscript data-raw/build_gene_files.R 2>&1 | tee build.log
```

The master script orchestrates nine steps:

1. UniProt human reference (~5 MB, ~1 min)
2. InterPro domains (~5 GB download, ~30 min)
3. ClinVar (~500 MB, ~10 min)
4. gnomAD v4.1 (per-gene GraphQL queries; ~6 hours for ~20K genes — adjust `sleep_per_query`)
5. AlphaMissense (1.2 GB, ~30 min — MD5 verified automatically)
6. REVEL (~700 MB, ~20 min)
7. CADD v1.7 (~300 GB download, can take ~12 hours; ~3 hours to parse)
8. Merge into per-gene bundles with validation
9. Manifest, licenses, README

Total wallclock: 12–24 hours. Most of this is downloads. Resumable
in the sense that each builder caches its source file — if the
script dies mid-run, restart and it skips already-downloaded files.

You'll get warnings about genes that failed validation. These are
written to `zenodo_payload/failed_genes.txt`. Spot-check a few; if
the failures look systematic (e.g. all from one source), debug the
relevant `build_*.R` script before proceeding.

---

## Step 3 — Sanity-check the output

```bash
ls zenodo_payload/ | head
wc -l zenodo_payload/MANIFEST.tsv   # should be ~20000
du -sh zenodo_payload/              # should be ~8-12 GB
```

Spot-check the largest gene file and a small one:

```r
# In R
big <- readRDS("zenodo_payload/TTN.rds")
str(big, max.level = 2)             # all seven elements present?
nrow(big$alphamissense)              # ~640,000 (TTN has 34K residues × 19 alts)
small <- readRDS("zenodo_payload/PATL1.rds")
str(small, max.level = 2)
```

Now test the package against the local files (without a Zenodo
upload yet):

```r
# Temporarily redirect the package's cache at the local payload
Sys.setenv(MSAVARIANT_CACHE = "/path/to/zenodo_payload_cache")
# Copy a gene file into the cache dir manually to simulate "downloaded"
dir.create("zenodo_payload_cache/0.1.0/genes", recursive = TRUE)
file.copy("zenodo_payload/PATL1.rds",
          "zenodo_payload_cache/0.1.0/genes/PATL1.rds")

library(msaVariant)
b <- fetch_gene_data("PATL1")
str(b, max.level = 2)
```

If `fetch_gene_data()` returns a valid bundle, the file format is
correct.

---

## Step 4 — Tarball and upload to Zenodo

```bash
cd zenodo_payload
tar --use-compress-program="xz -9 -T0" \
    -cvf ../msaVariant_data_v0.1.0.tar.xz .
cd ..
ls -lh msaVariant_data_v0.1.0.tar.xz   # expect ~6-9 GB
```

1. Go to https://zenodo.org/deposit/new
2. Fill in metadata:
   - **Title**: "msaVariant data deposit v0.1.0 — per-gene annotation bundles for human (InterPro, ClinVar, gnomAD v4.1, AlphaMissense, REVEL, CADD)"
   - **Authors**: Kaplan Lab members responsible for the build
   - **Description**: paste `zenodo_payload/README.md`
   - **Resource type**: Dataset
   - **License**: "Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC-BY-NC-SA-4.0)" *(required because of AlphaMissense + CADD)*
3. Upload `msaVariant_data_v0.1.0.tar.xz`
4. Click **Reserve DOI** (left sidebar). Zenodo reveals a DOI and
   a numeric record ID.
5. Click **Save** to keep your draft. **Don't publish yet** — go to
   step 5 first.

---

## Step 5 — Wire the DOI into the package

Edit `R/fetch_core.R`. Near the top:

```r
MSAVARIANT_DATA_VERSION <- "0.1.0"
MSAVARIANT_ZENODO_RECORD <- "1234567"           # the numeric ID
MSAVARIANT_DATA_DOI      <- "10.5281/zenodo.1234567"
```

Rebuild and test:

```bash
R CMD build .
R CMD check msaVariant_0.1.0.tar.gz
```

The check should pass. Manually try one fetch on a clean cache:

```r
# Clear any local cache
clear_cache()
b <- fetch_gene_data("PATL1")
stopifnot(!is.null(b))
```

This downloads from your draft Zenodo deposit. If it works, the
deposit is correctly configured.

---

## Step 6 — Publish (one-way!)

Go back to your Zenodo draft and click **Publish**. The record
becomes immutable and the DOI permanent.

Tag the package v0.1.0 on GitHub and push.

---

## Refreshing the data later

Quarterly or as needed:

1. `Rscript data-raw/build_gene_files.R` (rebuilds everything).
2. In Zenodo, open the concept record → "New version".
3. Upload the new tarball, reserve a new DOI.
4. Update `R/fetch_core.R` with the new DOI and bump
   `MSAVARIANT_DATA_VERSION`.
5. Bump the package version (e.g. v0.1.1) and publish a new Zenodo
   version + GitHub release.

Old DOIs stay live forever, so reproducibility of past analyses is
preserved.

---

## Troubleshooting

- **`build_clinvar.R` parse drops many rows**: HGVS strings in
  ClinVar are inconsistent. The current regex catches the common
  cases; ~5% are dropped silently. Acceptable for v0.1; revisit
  if a user reports a clinically important missing variant.

- **gnomAD API returns 429 (rate limited)**: increase
  `sleep_per_query` in `build_gnomad()` and rerun. The function
  is incremental — already-saved genes are skipped.

- **AlphaMissense MD5 mismatch**: the file was corrupted in
  transit. Delete `build_tmp/AlphaMissense_aa_substitutions.tsv.gz`
  and rerun.

- **CADD download dies**: use the European mirror by passing
  `european_mirror = TRUE` to `build_cadd()`. Resumable downloads
  via `wget -c` are recommended for partial-resume.

- **Validation warnings for individual genes**: the merger writes
  failed-gene names to `zenodo_payload/failed_genes.txt`. Inspect
  the warnings in `build.log` for the root cause. Common: a gene
  with an AlphaMissense UniProt accession that doesn't match our
  canonical (different isoform). Safe to skip these.
