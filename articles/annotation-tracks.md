# Annotation Tracks

Underneath the alignment,
[`plot_variant_overlay()`](https://ramizkrdnz.github.io/msaVariant/reference/plot_variant_overlay.md)
stacks up to six annotation tracks. Each one is a horizontal strip on
the same x coordinate system as the MSA, so a tile always sits directly
beneath the alignment column it describes. This article covers what each
track shows, where its data comes from, and — the part that is easy to
get wrong when reading a figure — how multiple records at one residue
are collapsed into a single tile.

``` r

library(msaVariant)

Sys.setenv(MSAVARIANT_CACHE = tempfile("msaVariant_cache_"))
import_local_bundle(
  system.file("extdata", "DEMO1.rds", package = "msaVariant"),
  gene = "DEMO1"
)
demo_fasta <- system.file("extdata", "demo_aligned.fasta",
                          package = "msaVariant")
demo <- fetch_gene_data("DEMO1")
```

> DEMO1 is synthetic demonstration data — fabricated tables in the real
> schema, not real ClinVar, gnomAD, AlphaMissense, REVEL or CADD
> records. Everything below is about *structure*, not about the biology
> of any real gene.

## The bundle

All six tracks read from one object: a **gene bundle**, a named list of
data frames fetched by
[`fetch_gene_data()`](https://ramizkrdnz.github.io/msaVariant/reference/fetch_gene_data.md)
or loaded from disk with
[`import_local_bundle()`](https://ramizkrdnz.github.io/msaVariant/reference/import_local_bundle.md).
There is no per-track network call.

``` r

names(demo)
#> [1] "meta"          "domains"       "clinvar"       "gnomad"       
#> [5] "alphamissense" "revel"         "cadd"
demo$meta[, c("gene", "uniprot_id", "protein_length")]
#>    gene uniprot_id protein_length
#> 1 DEMO1  DEMO00001             40
```

Each table also has an accessor. Note these take a **gene symbol**, not
a bundle object — they resolve the bundle from the cache themselves:

``` r

head(get_clinvar("DEMO1"), 3)
#>   pos aa_ref aa_alt aa_change      significance review_status  clinvar_id
#> 1   1      M      A       M1A Likely_pathogenic        2_star DEMO0000001
#> 2   2      K      A       K2A             Other        3_star DEMO0000002
#> 3   3      T      A       T3A     Likely_benign        1_star DEMO0000003
head(get_domains("DEMO1"), 3)
#>   start end                   name  accession   source
#> 1     5  35 Demo functional domain PFDEMO0001     Pfam
#> 2     2  10 Demo N-terminal region IPRDEMO001 InterPro
```

The full set:
[`get_gene_meta()`](https://ramizkrdnz.github.io/msaVariant/reference/get_gene_meta.md),
[`get_clinvar()`](https://ramizkrdnz.github.io/msaVariant/reference/get_clinvar.md),
[`get_gnomad()`](https://ramizkrdnz.github.io/msaVariant/reference/get_gnomad.md),
[`get_cadd()`](https://ramizkrdnz.github.io/msaVariant/reference/get_cadd.md),
[`get_revel()`](https://ramizkrdnz.github.io/msaVariant/reference/get_revel.md),
[`get_alphamissense()`](https://ramizkrdnz.github.io/msaVariant/reference/get_alphamissense.md),
[`get_domains()`](https://ramizkrdnz.github.io/msaVariant/reference/get_domains.md).
Use them when you want one table; hold the bundle from
[`fetch_gene_data()`](https://ramizkrdnz.github.io/msaVariant/reference/fetch_gene_data.md)
when you want several, to avoid re-resolving it each time.

## The reference figure

``` r

plot_variant_overlay(
  gene          = "DEMO1",
  aligned_fasta = demo_fasta,
  variant_pos   = 21,
  variant_label = "p.R21H"
)
```

![DEMO1 p.R21H overlay showing all six annotation
tracks](annotation-tracks_files/figure-html/figure-1.png)

Tracks appear top to bottom in a fixed order — domains, ClinVar, gnomAD,
CADD, AlphaMissense, REVEL — with the residue-position axis on the
lowest one and a legend dashboard at the foot.

## Domains

**Shows:** annotated functional regions as labelled boxes spanning their
residue range.

**Source:** InterPro and Pfam, carried in `bundle$domains` as `start`,
`end`, `name`, `accession`, `source`.

``` r

get_domains("DEMO1")
#>   start end                   name  accession   source
#> 1     5  35 Demo functional domain PFDEMO0001     Pfam
#> 2     2  10 Demo N-terminal region IPRDEMO001 InterPro
```

**Aggregation:** none — domains are ranges, not per-residue values, and
are drawn as rectangles from `start` to `end`. What *is* selective is
the **label**. Overlapping domains are common, and stacking every name
in a one-line strip is unreadable, so one domain is chosen to label:
Pfam entries are preferred over generic superfamily-style names, then a
domain containing the variant wins over one that merely overlaps the
window. The boxes for the other visible domains are still drawn.

Here, residue 21 falls inside the 5–35 Pfam “Demo functional domain”,
which is what gets the label — and what fires PM1.

## ClinVar

**Shows:** clinical significance assertions, one coloured tile per
annotated residue.

**Source:** ClinVar, in `bundle$clinvar` as `pos`, `aa_ref`, `aa_alt`,
`aa_change`, `significance`, `review_status`, `clinvar_id`.

``` r

head(get_clinvar("DEMO1"), 5)
#>   pos aa_ref aa_alt aa_change      significance review_status  clinvar_id
#> 1   1      M      A       M1A Likely_pathogenic        2_star DEMO0000001
#> 2   2      K      A       K2A             Other        3_star DEMO0000002
#> 3   3      T      A       T3A     Likely_benign        1_star DEMO0000003
#> 4   4      A      G       A4G       Conflicting        3_star DEMO0000004
#> 5   5      Y      A       Y5A            Benign        1_star DEMO0000005
table(get_clinvar("DEMO1")$significance)
#> 
#>        Pathogenic Likely_pathogenic               VUS     Likely_benign 
#>                 4                 8                 2                 4 
#>            Benign       Conflicting             Other 
#>                 5                 4                 2
```

**Aggregation:** this is the one categorical track, and it is *not*
numerically aggregated — each row is drawn as a tile at its position,
coloured by `significance` through a named palette covering
`Pathogenic`, `Likely_pathogenic`, `VUS`, `Likely_benign`, `Benign`,
`Conflicting` and `Other`. Where several records share a residue, they
are drawn at the same column and the last one plotted is what you see on
top.

So read the ClinVar strip as *“something is reported here”*, and go to
the table when the exact assertion matters:

``` r

subset(get_clinvar("DEMO1"), pos == 21)
#>    pos aa_ref aa_alt aa_change      significance review_status  clinvar_id
#> 15  21      R      A      R21A        Pathogenic        4_star DEMO0000015
#> 28  21      R      H      R21H        Pathogenic        2_star DEMO9000001
#> 29  21      R      C      R21C Likely_pathogenic        1_star DEMO9000002
```

Three records at residue 21 — which is exactly what lets both PS1 (same
change, pathogenic) and PM5 (different change, same residue, pathogenic)
fire. See
[`vignette("acmg-evidence-codes")`](https://ramizkrdnz.github.io/msaVariant/articles/acmg-evidence-codes.md).

## gnomAD

**Shows:** population allele frequency, as a sequential colour ramp.

**Source:** gnomAD, in `bundle$gnomad` as `pos`, `aa_change`,
`consequence`, `af_joint`, `ac_joint`, `an_joint`, `filter`.

``` r

head(get_gnomad("DEMO1"), 5)
#>   pos aa_ref aa_alt aa_change consequence  af_joint ac_joint an_joint filter
#> 1   1      M      A       M1A    missense 0.0005818       12   152312   PASS
#> 2   5      Y      A       Y5A  synonymous 0.0009313       26   152312   PASS
#> 3   6      I      A       I6A    missense 0.0036489        1   152312   PASS
#> 4   8      K      A       K8A  synonymous 0.0020599        2   152312   PASS
#> 5   9      Q      A       Q9A    missense 0.0020708       24   152312   PASS
```

**Aggregation:** the **maximum** `af_joint` per position, then
[`log10()`](https://rdrr.io/r/base/Log.html). Maximum rather than mean
because the question the track answers is “how common is variation at
this residue at all?” — one common substitution is the
constraint-relevant signal, and averaging it against rare neighbours
would hide it. The log transform is what makes the ramp readable: allele
frequencies span orders of magnitude, and on a linear scale every rare
variant collapses into one indistinguishable colour.

Residues absent from gnomAD get no tile — a gap in the strip means “not
observed”, which is itself the PM2 evidence.

``` r

nrow(subset(get_gnomad("DEMO1"), pos == 21))
#> [1] 0
```

## CADD

**Shows:** CADD PHRED-scaled deleteriousness, on a fixed 0–40 scale.

**Source:** CADD, in `bundle$cadd` as `pos`, `aa_change`, `consequence`,
`cadd_raw`, `cadd_phred`.

``` r

head(get_cadd("DEMO1"), 5)
#>   pos aa_ref aa_alt aa_change consequence cadd_raw cadd_phred
#> 1   1      M      A       M1A    missense    1.535      12.28
#> 2   2      K      A       K2A    missense    0.539       4.31
#> 3   3      T      A       T3A    missense    4.897      39.17
#> 4   4      A      G       A4G    missense    2.485      19.88
#> 5   5      Y      A       Y5A    missense    0.465       3.72
```

**Aggregation:** the **median** `cadd_phred` per position. A residue can
carry several possible substitutions with quite different scores; the
median gives the typical consequence of mutating that residue without
letting a single extreme substitution dominate the tile.

The scale is pinned to `c(0, 40)` rather than fitted to the data, so
colour means the same thing across every figure you make. PHRED 20 is
the top 1% of scored substitutions, 30 the top 0.1%; the PP3 threshold
in this package is 25.

## AlphaMissense

**Shows:** AlphaMissense pathogenicity score, fixed 0–1 scale.

**Source:** AlphaMissense, in `bundle$alphamissense` as `pos`,
`aa_change`, `am_score`, `am_class`.

``` r

head(get_alphamissense("DEMO1"), 5)
#>   pos aa_ref aa_alt aa_change am_score      am_class
#> 1   1      M      A       M1A   0.2457 likely_benign
#> 2   2      K      A       K2A   0.3511     ambiguous
#> 3   3      T      A       T3A   0.1590 likely_benign
#> 4   4      A      G       A4G   0.3041 likely_benign
#> 5   5      Y      A       Y5A   0.0175 likely_benign
table(get_alphamissense("DEMO1")$am_class)
#> 
#>     likely_benign         ambiguous likely_pathogenic 
#>                19                10                11
```

**Aggregation:** the **mean** `am_score` per position. Mean rather than
median here because AlphaMissense scores are bounded, roughly
continuous, and typically dense across substitutions at a residue, so
the average reads as a per-residue tolerance estimate.

Note the track draws the continuous `am_score`, while the ACMG PP3 rule
reads the categorical `am_class`. A residue can look mid-range in the
strip while a specific substitution at it is still classified
`likely_pathogenic` — the strip summarises the residue, the code
evaluates the substitution.

## REVEL

**Shows:** REVEL ensemble missense score, fixed 0–1 scale.

**Source:** REVEL, in `bundle$revel` as `pos`, `aa_change`,
`revel_score`.

``` r

head(get_revel("DEMO1"), 5)
#>   pos aa_ref aa_alt aa_change revel_score
#> 1   1      M      A       M1A      0.8868
#> 2   2      K      A       K2A      0.1363
#> 3   3      T      A       T3A      0.7853
#> 4   4      A      G       A4G      0.4533
#> 5   5      Y      A       Y5A      0.1357
```

**Aggregation:** the **mean** `revel_score` per position, same reasoning
as AlphaMissense. The PP3 threshold is `> 0.7`.

REVEL and AlphaMissense are both missense-pathogenicity predictors and
often agree, which is deliberate: the point of showing both is that
their *disagreements* are informative, and `pp3_min_predictors` lets you
decide how much agreement you require.

## Summary

| Track         | Data            | Per-residue value drawn   | Scale       |
|---------------|-----------------|---------------------------|-------------|
| Domains       | InterPro / Pfam | range boxes, one labelled | —           |
| ClinVar       | ClinVar         | `significance` per record | categorical |
| gnomAD        | gnomAD          | `log10(max(af_joint))`    | fitted      |
| CADD          | CADD            | `median(cadd_phred)`      | fixed 0–40  |
| AlphaMissense | AlphaMissense   | `mean(am_score)`          | fixed 0–1   |
| REVEL         | REVEL           | `mean(revel_score)`       | fixed 0–1   |

The recurring caveat: **tracks are per-residue, ACMG codes are
per-substitution.** The strips are for spotting where in the alignment
the evidence clusters;
[`compute_acmg_codes()`](https://ramizkrdnz.github.io/msaVariant/reference/compute_acmg_codes.md)
is what evaluates the specific change.

## Coordinate mapping

Bundle tables are in reference protein (UniProt canonical) numbering.
Alignments have gaps, so column 30 of the MSA is generally not residue
30.
[`map_variant_to_msa()`](https://ramizkrdnz.github.io/msaVariant/reference/map_variant_to_msa.md)
does the translation, and every track is mapped through it before
drawing:

``` r

map_variant_to_msa(c(5, 21, 35), demo_fasta, ref_name = "DEMO1_HUMAN")
#> [1]  5 21 35
```

[`build_msa_coord_map()`](https://ramizkrdnz.github.io/msaVariant/reference/build_msa_coord_map.md)
returns the full residue-to-column table if you need to place your own
annotations by hand.

Positions that fall on a gap in the reference sequence map to `NA` and
are dropped rather than drawn at the wrong column.

## Bringing your own

Any per-residue value you can put in a data frame with a `pos` column
can be added as a track with
[`geom_track()`](https://ramizkrdnz.github.io/msaVariant/reference/geom_track.md),
and a bundle you build yourself can be loaded with
[`import_local_bundle()`](https://ramizkrdnz.github.io/msaVariant/reference/import_local_bundle.md).
See
[`vignette("bring-your-own-data")`](https://ramizkrdnz.github.io/msaVariant/articles/bring-your-own-data.md).

## See also

- [`vignette("acmg-evidence-codes")`](https://ramizkrdnz.github.io/msaVariant/articles/acmg-evidence-codes.md)
  — the rules these tables feed
- [`vignette("toggling-layers")`](https://ramizkrdnz.github.io/msaVariant/articles/toggling-layers.md)
  — hiding tracks you do not need
- [`vignette("colour-schemes")`](https://ramizkrdnz.github.io/msaVariant/articles/colour-schemes.md)
  — restyling track palettes
