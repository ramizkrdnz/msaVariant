# msaVariant 0.2.0

## New features

* `plot_variant_overlay()` assembles a complete publication-ready figure:
  a cross-species MSA window, per-residue annotation tracks (ClinVar,
  gnomAD, AlphaMissense, REVEL, CADD, protein domains), and an automatic
  ACMG evidence-code strip (PS1, PM1, PM2, PM5, PP3).
* `compute_acmg_codes()` derives the ACMG evidence codes for a single
  substitution entirely from a gene bundle, with a user-tunable
  `pp3_min_predictors` threshold.
* Per-figure customisation: three colour schemes (`journal`,
  `colorblind`, `grayscale`) plus per-element colour overrides, and an
  independent on/off toggle for every annotation layer with automatic
  layout re-flow.
* Data layer: per-gene annotation bundles fetched from a Zenodo deposit
  (`fetch_gene_data()`), cached locally, with a `MANIFEST.tsv` integrity
  check and `available_genes()`. Bundles can also be built locally and
  loaded with `import_local_bundle()`.
* `geom_track()` overlays arbitrary per-residue data on an alignment for
  bring-your-own-data workflows.
* Ships a synthetic, licence-clean `DEMO1` example bundle so the
  documentation and examples run out of the box.

## Documentation

* Seven vignettes: a Get Started guide, a full workflow tutorial, and
  guides for ACMG evidence codes, annotation tracks, colour schemes,
  layer toggling, and bring-your-own-data.
