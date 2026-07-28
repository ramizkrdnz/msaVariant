# msaVariant: Clinical-genetics MSA visualisation with variant overlay

Extends \`ggmsa\` with annotation layers tailored to clinical and
rare-disease genetics workflows.

## Details

What the package provides (v0.1):

\* \`geom_variant()\` — user variants (tidy data frame) \*
\`geom_domain()\` — protein domain track (user-supplied) \*
\`geom_track()\` — generic per-residue annotation track for any
user-supplied data (gnomAD AF, ClinVar significance, AlphaMissense
scores, PTM sites, etc.) \* \`conservation_score()\` — per-column
Shannon / Jensen-Shannon computed from the MSA itself

Plus utility functions:

\* \`build_msa_coord_map()\` — residue-to-column lookup \*
\`map_variant_to_msa()\` — vectorised position mapping \*
\`scale_fill_pathogenicity()\` — colourblind-safe ACMG palette

What the package does NOT bundle:

Population-genetics and pathogenicity databases (gnomAD, ClinVar,
AlphaMissense, REVEL, CADD) change frequently, are licensed in ways that
constrain redistribution, and are too large to ship as static files.
Users supply their own annotation data via \`geom_track()\`. A future
v0.2 will add optional online fetchers that pull these databases on
demand and cache results locally.

## See also

Useful links:

- <https://github.com/thekaplanlab/msaVariant>

- <https://thekaplanlab.github.io/msaVariant/>

- Report bugs at <https://github.com/thekaplanlab/msaVariant/issues>

## Author

**Maintainer**: Ramiz Karadeniz <ramiz.karadeniz81@gmail.com>
([ORCID](https://orcid.org/0009-0000-3290-4707))

Authors:

- Ramiz Karadeniz <ramiz.karadeniz81@gmail.com>
  ([ORCID](https://orcid.org/0009-0000-3290-4707))

- Oktay I. Kaplan ([ORCID](https://orcid.org/0000-0002-8733-0920))

- Sebiha Cevik ([ORCID](https://orcid.org/0000-0002-0935-1929))

Other contributors:

- Guangchuang Yu (ggmsa, on which this package builds) \[contributor\]
