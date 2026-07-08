% Generated manually to accompany R/compute_acmg_codes.R
\name{compute_acmg_codes}
\alias{compute_acmg_codes}
\title{Compute triggered ACMG evidence codes for one variant}
\usage{
compute_acmg_codes(bundle, variant_pos, aa_change, pp3_min_predictors = 2L)
}
\arguments{
\item{bundle}{A gene bundle: a named list with elements \code{domains},
\code{clinvar}, \code{gnomad}, \code{alphamissense}, \code{revel},
\code{cadd} (as produced by \code{\link{import_local_bundle}} /
\code{\link{fetch_gene_data}}).}

\item{variant_pos}{Integer residue position (canonical UniProt numbering).}

\item{aa_change}{Character variant identifier, e.g. \code{"R175H"} or
\code{"p.R175H"}. Ref/alt are parsed from this string.}

\item{pp3_min_predictors}{Integer 1-3 (default 2). PP3 fires when at least
this many of the three PP3 predictors (AlphaMissense class ==
"likely_pathogenic", REVEL > 0.7, CADD phred >= 25) pass their threshold.}
}
\value{
A character vector of the triggered codes, in canonical
evidence-strength order (\code{PS1}, \code{PM1}, \code{PM2}, \code{PM5},
\code{PP3}); \code{character(0)} if none fire.
}
\description{
Evaluates a fixed subset of ACMG/AMP evidence codes for a single
missense substitution using only the annotation tables carried in a gene
bundle. Returns just the codes that are triggered.
}
\details{
Rules (all computed from the bundle; no external lookups):
\describe{
  \item{PS1}{ClinVar contains the exact same \code{aa_change} with
    significance \code{"Pathogenic"}.}
  \item{PM1}{\code{variant_pos} falls within the start-end range of at
    least one \code{domains} row.}
  \item{PM2}{This \code{aa_change} is absent from gnomAD, or its
    \code{af_joint} is below 0.0001.}
  \item{PM5}{ClinVar has a record at the same position with a different
    alt residue and significance \code{"Pathogenic"} or
    \code{"Likely_pathogenic"}.}
  \item{PP3}{Count how many of the three predictors pass — AlphaMissense
    \code{am_class == "likely_pathogenic"}, REVEL \code{revel_score > 0.7},
    CADD \code{cadd_phred >= 25}; triggered when the count is at least
    \code{pp3_min_predictors}.}
}
}
\examples{
\dontrun{
b <- fetch_gene_data("TP53")
compute_acmg_codes(b, 175, "R175H")
}
}
