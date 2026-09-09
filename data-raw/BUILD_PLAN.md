# msaVariant — real grouped-deposit BUILD PLAN

Concrete plan for producing the real ~19,500-gene grouped data deposit,
with a **100–500-gene pilot first**. Decisions locked in by the maintainer:

- **CADD + REVEL come from Ensembl REST per-gene** (~7 GB total footprint),
  NOT the 300 GB bulk CADD file. This mirrors the provenance that produced
  the proven TP53/WDR31 bundles.
- **Pilot before full run.** The automated pipeline has never produced a
  real bundle end-to-end; validate at small scale first.

---

## 0. Why a pilot — the critical context

The real TP53/WDR31 bundles were built by `build_gene_bundle.R` (repo root)
from **manually-downloaded per-gene web exports** — that path is proven but
does not scale. The automated pipeline (`data-raw/build_scripts/*` +
`build_gene_files.R`) is the scalable design but was **never run to produce a
real bundle**. The pilot de-risks it before any multi-hour crawl.

---

## 1. Source strategy & footprint

| Source | Fetch mode | Footprint | Notes |
|---|---|---|---|
| UniProt ref | full REST TSV (reviewed human) | ~5 MB | gene list origin (§2) |
| ClinVar | full `variant_summary.txt.gz` | ~500 MB | cached; slice by gene |
| InterPro/domains | full `protein2ipr.dat.gz` | ~5 GB | cached; streamed |
| AlphaMissense | full `.tsv.gz` (MD5-checked) | 1.2 GB | cached; streamed |
| gnomAD | **per-gene GraphQL API** | — (~6–8 h) | now checkpointed (§3) |
| **CADD + REVEL** | **per-gene Ensembl REST** | streamed slices | **new builder (§4)** — replaces bulk CADD (300 GB) + REVEL (700 MB) |

**Realistic minimum: ~7 GB download, ~15 GB disk** (vs ~300 GB / ~500 GB with
bulk CADD). CADD/REVEL and gnomAD are per-gene; only ClinVar/InterPro/
AlphaMissense/UniProt are bulk-downloaded once and cached.

---

## 2. Gene list

Derived from **UniProt** (`build_uniprot_ref.R`): reviewed SwissProt human
(`organism_id:9606 AND reviewed:true`), kept where a gene symbol and valid
protein length exist → **~19,500** canonical entries (isoforms stripped). This
is UniProt-reviewed protein-coding, not the full HGNC set. The merge iterates
exactly this reference.

---

## 3. Resumability — DONE (this change)

A 19,500-gene, 6–8 h run cannot restart from zero. Hardened:

- **`merge_into_gene_bundles.R`**: per-gene `tryCatch` (a hard error on one
  gene logs + continues, never halts); `if (file.exists(<gene>.rds)) next`
  resume skip-guard; **atomic write** (temp + rename) so a crash never leaves
  a truncated file the skip-guard would mistake for done. Returns a
  written/skipped/failed summary; failures → `failed_genes.txt`.
- **`build_gnomad.R`**: genuine per-gene checkpoint cache
  (`<tmp>/gnomad_cache/<gene>.rds`). A definitive 200 response (even 0 rows)
  is cached so it is never re-queried; **transient failures are NOT cached**
  so they retry next run. Per-gene `tryCatch`. This makes the code match what
  `ZENODO_UPLOAD.md` already claimed ("incremental — already-saved genes are
  skipped") — previously the doc was aspirational; now it is true.
- The **Ensembl builder (§4)** uses the same checkpoint + tryCatch pattern.

Known pre-existing bug to fix during the pilot: `build_gene_files.R` calls
`build_manifest(out_dir=...)` but the parameter is `payload_dir=`. Use
`payload_dir=`.

---

## 4. Ensembl CADD+REVEL builder — DRAFTED (needs pilot validation)

`data-raw/build_scripts/build_ensembl_cadd_revel.R` — one VEP query per gene
yields both CADD and REVEL; returns `list(cadd = <by-gene>, revel = <by-gene>)`
and **replaces** `build_cadd()` + `build_revel()` in the orchestrator. Parsing
semantics are ported verbatim from the proven manual builder
(`build_tp53_from_real_data.R` §5): AA `"K/R"` → ref/alt, protein position →
`pos`, max PHRED / max REVEL per substitution, REVEL missense-only, clip to
protein length. Resumable + tryCatch.

**⚠️ The pilot's FIRST job:** confirm the public Ensembl REST actually returns
**both** scores. CADD is available via `?CADD=1`; REVEL is served via the
dbNSFP plugin (`?dbNSFP=REVEL_score`) which may **not** be enabled on the
public server, and the exact VEP-by-gene endpoint + response field names
(`transcript_consequences`, `cadd_phred`, `cadd_raw`, `revel_score`,
`amino_acids`, `protein_start`) must be checked against a live TP53 response.
**Validation gate:** run the builder on TP53 alone and diff its cadd/revel
tables against the known-good `tests/testthat/fixtures/TP53.rds`. If REVEL is
absent from REST, fall back (in order): (a) dbNSFP tabix range-slice per gene,
(b) REVEL bulk file + CADD from REST. **Do not launch the full crawl until
this passes.**

---

## 5. Grouping — DONE (this change)

`data-raw/build_group_bundles.R` — `build_group_bundles(gene_dir, out_dir, ...)`
consumes the merge output and packs per-gene bundles into groups with a **dual
cap**: a new group starts when adding the next gene would exceed **either**
`max_genes` (~200) **or** `max_bytes` (~40 MB). Alphabetical greedy →
deterministic; a single gene over the byte cap (e.g. TTN) gets its own group.
Writes `group_NNN.rds` (named lists of bundles), auto-generates
`gene_group_index.tsv` (shipped to `inst/extdata/` + a self-describing copy in
the payload), and runs `build_manifest()`. Validated end-to-end: real fixtures
pack correctly (count- and byte-caps both fire), and the output round-trips
through the runtime `fetch_gene_data()` + `available_genes()`.

Runtime needs **no change** — the grouped fetcher, auto-manifest download, and
`available_genes()` already handle any scale.

---

## 6. Upload to Zenodo

~98 group files + MANIFEST. Use **`zen4R`** (R, matches the pipeline; install
it). Rehearse on **sandbox** first (`zen4R` `sandbox = TRUE`), exactly as the
prototype was validated. Flow: create deposition → upload group files +
MANIFEST → metadata (**CC-BY-NC-SA-4.0**, required by AlphaMissense/CADD) →
**Reserve DOI** → publish. Then wire `MSAVARIANT_ZENODO_RECORD` (numeric id) +
`MSAVARIANT_DATA_DOI` in `R/fetch_core.R` and bump the data/package version.

---

## 7. Pilot procedure (next, ~100 genes incl. BRCA1/TTN)

1. `build_uniprot_ref()` → take ~100 genes, **including large ones (BRCA1,
   TTN)** to exercise size-balancing, plus TP53/WDR31 as ground truth.
2. **Validate the Ensembl builder on TP53** vs the fixture (§4 gate) before
   crawling the rest.
3. Run the 6 source builders (CADD/REVEL via Ensembl) → intermediates.
4. `merge_into_gene_bundles()` → per-gene `.rds` (resume-safe).
5. `build_group_bundles()` → group files + index + MANIFEST.
6. Upload the pilot groups to a **sandbox** record; run the grouped E2E
   (download-once, cache-hit, checksum, `available_genes()`, render) exactly
   as done for record 596247.
7. Only after the pilot is green: full ~19,500 run, then production Zenodo.

**Checkpoint to review before the full run:** pilot bundle spot-checks
(TP53/BRCA1 match expectations), group size distribution, `failed_genes.txt`
contents, and total wall-clock extrapolation for gnomAD + Ensembl at 19,500.
