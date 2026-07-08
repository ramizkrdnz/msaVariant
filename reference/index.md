# Package index

## Plotting

Assemble a complete publication-ready variant overlay figure.

- [`plot_variant_overlay()`](https://ramizkrdnz.github.io/msaVariant/reference/plot_variant_overlay.md)
  : Generate a publication-grade variant-overlay figure

## Annotation geoms

ggplot2 layers that paint variants and per-residue annotation tracks
onto a multiple sequence alignment.

- [`geom_variant()`](https://ramizkrdnz.github.io/msaVariant/reference/geom_variant.md)
  : Overlay variants onto a multiple sequence alignment
- [`geom_track()`](https://ramizkrdnz.github.io/msaVariant/reference/geom_track.md)
  : Overlay an arbitrary per-residue annotation track onto an MSA
- [`geom_domain()`](https://ramizkrdnz.github.io/msaVariant/reference/geom_domain.md)
  : Overlay protein domain annotation on an MSA
- [`geom_clinvar()`](https://ramizkrdnz.github.io/msaVariant/reference/geom_clinvar.md)
  : Overlay ClinVar pathogenicity calls onto an MSA
- [`geom_gnomad()`](https://ramizkrdnz.github.io/msaVariant/reference/geom_gnomad.md)
  : Overlay gnomAD allele frequencies onto an MSA
- [`geom_cadd()`](https://ramizkrdnz.github.io/msaVariant/reference/geom_cadd.md)
  : Overlay CADD per-residue PHRED scores onto an MSA
- [`geom_alphamissense()`](https://ramizkrdnz.github.io/msaVariant/reference/geom_alphamissense.md)
  : Overlay AlphaMissense per-residue mean scores onto an MSA
- [`geom_revel()`](https://ramizkrdnz.github.io/msaVariant/reference/geom_revel.md)
  : Overlay REVEL per-residue scores onto an MSA

## ACMG evidence codes

Automatic ACMG/AMP variant-classification evidence codes.

- [`compute_acmg_codes()`](https://ramizkrdnz.github.io/msaVariant/reference/compute_acmg_codes.md)
  : Compute triggered ACMG evidence codes for one variant

## Colours & scales

Publication-quality colour schemes and pathogenicity fills.

- [`resolve_plot_colors()`](https://ramizkrdnz.github.io/msaVariant/reference/resolve_plot_colors.md)
  : Resolve the full color specification for a variant-overlay figure
- [`scale_fill_pathogenicity()`](https://ramizkrdnz.github.io/msaVariant/reference/scale_fill_pathogenicity.md)
  : Default colour scale for ACMG pathogenicity tiers

## Conservation & coordinate mapping

Per-column conservation scoring and variant-to-MSA coordinate mapping.

- [`conservation_score()`](https://ramizkrdnz.github.io/msaVariant/reference/conservation_score.md)
  : Per-column conservation score for an MSA
- [`map_variant_to_msa()`](https://ramizkrdnz.github.io/msaVariant/reference/map_variant_to_msa.md)
  : Map a vector of variant residue positions to MSA columns
- [`build_msa_coord_map()`](https://ramizkrdnz.github.io/msaVariant/reference/build_msa_coord_map.md)
  : Build a residue -\> MSA column lookup

## Data fetching & cache

Retrieve per-gene annotation bundles and manage the local cache.

- [`fetch_gene_data()`](https://ramizkrdnz.github.io/msaVariant/reference/fetch_gene_data.md)
  : Fetch the combined annotation file for a gene
- [`available_genes()`](https://ramizkrdnz.github.io/msaVariant/reference/available_genes.md)
  : List genes available in the local data manifest
- [`clear_cache()`](https://ramizkrdnz.github.io/msaVariant/reference/clear_cache.md)
  : Clear the msaVariant download cache
- [`cache_location()`](https://ramizkrdnz.github.io/msaVariant/reference/cache_location.md)
  : Show the location of the msaVariant cache
- [`cache_summary()`](https://ramizkrdnz.github.io/msaVariant/reference/cache_summary.md)
  : Summarize what's currently cached

## Data accessors

Pull individual annotation tables out of a fetched gene bundle.

- [`get_gene_meta()`](https://ramizkrdnz.github.io/msaVariant/reference/get_gene_meta.md)
  : Get gene metadata
- [`get_clinvar()`](https://ramizkrdnz.github.io/msaVariant/reference/get_clinvar.md)
  : Get ClinVar variants for a gene
- [`get_gnomad()`](https://ramizkrdnz.github.io/msaVariant/reference/get_gnomad.md)
  : Get gnomAD per-variant allele frequencies for a gene
- [`get_cadd()`](https://ramizkrdnz.github.io/msaVariant/reference/get_cadd.md)
  : Get CADD per-substitution PHRED scores for a gene
- [`get_revel()`](https://ramizkrdnz.github.io/msaVariant/reference/get_revel.md)
  : Get REVEL per-substitution scores for a gene
- [`get_alphamissense()`](https://ramizkrdnz.github.io/msaVariant/reference/get_alphamissense.md)
  : Get AlphaMissense per-substitution scores for a gene
- [`get_domains()`](https://ramizkrdnz.github.io/msaVariant/reference/get_domains.md)
  : Get InterPro/Pfam domains for a gene

## Data structures & validation

Construct and validate the gene-bundle data structure.

- [`validate_gene_data()`](https://ramizkrdnz.github.io/msaVariant/reference/validate_gene_data.md)
  : Validate a per-gene annotation object against the package spec
- [`empty_gene_data()`](https://ramizkrdnz.github.io/msaVariant/reference/empty_gene_data.md)
  : Construct an empty per-gene object conforming to the spec

## Import

Load a locally built gene bundle without a remote fetch.

- [`import_local_bundle()`](https://ramizkrdnz.github.io/msaVariant/reference/import_local_bundle.md)
  : Import a locally-built gene bundle into the cache

## Package overview

- [`msaVariant-package`](https://ramizkrdnz.github.io/msaVariant/reference/msaVariant-package.md)
  [`msaVariant`](https://ramizkrdnz.github.io/msaVariant/reference/msaVariant-package.md)
  : msaVariant: Clinical-genetics MSA visualisation with variant overlay
