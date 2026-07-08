## -----------------------------------------------------------------
## msaVariant: per-gene data schema validator
## -----------------------------------------------------------------
##
## Enforces the contract defined in DATA_FORMAT_SPEC.md.
## Called automatically by `fetch_gene_data()` after each download
## and after each cache hit. Also usable directly to validate a
## locally-built file before uploading to Zenodo.

# Canonical factor levels per the spec
.SIGNIFICANCE_LEVELS <- c("Pathogenic", "Likely_pathogenic", "VUS",
                          "Likely_benign", "Benign",
                          "Conflicting", "Other")
.REVIEW_STATUS_LEVELS <- c("4_star", "3_star", "2_star", "1_star", "0_star")
.DOMAIN_SOURCE_LEVELS <- c("Pfam", "InterPro", "SMART", "PROSITE",
                           "PRINTS", "PANTHER")
.CONSEQUENCE_LEVELS <- c("missense", "synonymous", "stop_gained",
                         "frameshift", "inframe_deletion",
                         "inframe_insertion", "splice_donor",
                         "splice_acceptor", "other")
.AM_CLASS_LEVELS <- c("likely_benign", "ambiguous", "likely_pathogenic")

# Required column schemas (NULL type = "any")
.SCHEMA_META <- list(
  gene                  = "character",
  uniprot_id            = "character",
  protein_length        = "integer",
  ensembl_gene_id       = "character",
  ensembl_transcript_id = "character",
  build_date            = "Date",
  source_versions       = "character"
)
.SCHEMA_DOMAINS <- list(
  start     = "integer",
  end       = "integer",
  name      = "character",
  accession = "character",
  source    = "factor"
)
.SCHEMA_CLINVAR <- list(
  pos            = "integer",
  aa_ref         = "character",
  aa_alt         = "character",
  aa_change      = "character",
  significance   = "factor",
  review_status  = "factor",
  clinvar_id     = "character"
)
.SCHEMA_GNOMAD <- list(
  pos         = "integer",
  aa_ref      = "character",
  aa_alt      = "character",
  aa_change   = "character",
  consequence = "factor",
  af_joint    = "numeric",
  ac_joint    = "integer",
  an_joint    = "integer",
  filter      = "character"
)
.SCHEMA_ALPHAMISSENSE <- list(
  pos        = "integer",
  aa_ref     = "character",
  aa_alt     = "character",
  aa_change  = "character",
  am_score   = "numeric",
  am_class   = "factor"
)
.SCHEMA_REVEL <- list(
  pos         = "integer",
  aa_ref      = "character",
  aa_alt      = "character",
  aa_change   = "character",
  revel_score = "numeric"
)
.SCHEMA_CADD <- list(
  pos         = "integer",
  aa_ref      = "character",
  aa_alt      = "character",
  aa_change   = "character",
  consequence = "factor",
  cadd_raw    = "numeric",
  cadd_phred  = "numeric"
)

.REQUIRED_TABLES <- c("meta", "domains", "clinvar", "gnomad",
                      "alphamissense", "revel", "cadd")

.check_df_schema <- function(df, schema, table_name, factor_levels = list()) {
  errs <- character(0)
  if (!is.data.frame(df)) {
    errs <- c(errs, sprintf("%s: not a data.frame.", table_name))
    return(errs)
  }
  missing <- setdiff(names(schema), names(df))
  if (length(missing)) {
    errs <- c(errs, sprintf("%s: missing required columns: %s",
                            table_name, paste(missing, collapse = ", ")))
  }
  # Type checks only on columns that exist; missing ones already reported
  for (col in intersect(names(schema), names(df))) {
    expected <- schema[[col]]
    actual   <- if (is.factor(df[[col]])) "factor" else typeof(df[[col]])
    ok <- switch(
      expected,
      "integer"   = is.integer(df[[col]]) || is.numeric(df[[col]]),
      "numeric"   = is.numeric(df[[col]]),
      "character" = is.character(df[[col]]),
      "factor"    = is.factor(df[[col]]),
      "Date"      = inherits(df[[col]], "Date"),
      TRUE
    )
    if (!ok) {
      errs <- c(errs, sprintf("%s$%s: expected %s, got %s",
                              table_name, col, expected, actual))
    }
    # Factor-levels check: must match canonical set exactly
    if (expected == "factor" && is.factor(df[[col]]) &&
        !is.null(factor_levels[[col]])) {
      if (!identical(levels(df[[col]]), factor_levels[[col]])) {
        errs <- c(errs, sprintf("%s$%s: factor levels mismatch (got: %s)",
                                table_name, col,
                                paste(levels(df[[col]]), collapse = ",")))
      }
    }
  }
  # No NA in required columns
  for (col in intersect(names(schema), names(df))) {
    if (nrow(df) > 0L && anyNA(df[[col]])) {
      errs <- c(errs, sprintf("%s$%s: NA in required column",
                              table_name, col))
    }
  }
  errs
}

#' Validate a per-gene annotation object against the package spec
#'
#' Checks that the object conforms to the structure defined in
#' `DATA_FORMAT_SPEC.md`: presence of all 7 named elements, required
#' columns, correct types, and no `NA` in required columns. Returns
#' a list with `valid` (logical) and `message` (character) describing
#' all issues found.
#'
#' Used internally by `fetch_gene_data()` after every cache hit or
#' download. Also useful for the lab's build pipeline to validate
#' files before uploading to Zenodo.
#'
#' @param obj The deserialized object from a per-gene `.rds` file.
#' @param strict If `TRUE`, throw an error on validation failure
#'   instead of returning a list. Default `FALSE`.
#' @return If `strict = FALSE`: a list with elements `valid`
#'   (logical) and `message` (character vector of issues, or NULL
#'   if valid). If `strict = TRUE`: invisibly `TRUE` on success, or
#'   `stop()` with a formatted message on failure.
#' @export
validate_gene_data <- function(obj, strict = FALSE) {
  errs <- character(0)

  if (!is.list(obj) || is.data.frame(obj)) {
    return(list(valid = FALSE,
                message = "Object is not a list."))
  }

  # 1. All required tables present
  missing_tables <- setdiff(.REQUIRED_TABLES, names(obj))
  if (length(missing_tables)) {
    for (t in missing_tables) {
      errs <- c(errs, sprintf("Missing top-level element: %s", t))
    }
  }

  # 2. Per-table schema checks
  if ("meta" %in% names(obj)) {
    errs <- c(errs, .check_df_schema(obj$meta, .SCHEMA_META, "meta"))
    if (is.data.frame(obj$meta) && nrow(obj$meta) != 1L) {
      errs <- c(errs, sprintf("meta: must have exactly 1 row (has %d).",
                              nrow(obj$meta)))
    }
  }
  if ("domains" %in% names(obj)) {
    errs <- c(errs, .check_df_schema(obj$domains, .SCHEMA_DOMAINS, "domains",
                                     factor_levels = list(source = .DOMAIN_SOURCE_LEVELS)))
  }
  if ("clinvar" %in% names(obj)) {
    errs <- c(errs, .check_df_schema(obj$clinvar, .SCHEMA_CLINVAR, "clinvar",
                                     factor_levels = list(
                                       significance = .SIGNIFICANCE_LEVELS,
                                       review_status = .REVIEW_STATUS_LEVELS)))
  }
  if ("gnomad" %in% names(obj)) {
    errs <- c(errs, .check_df_schema(obj$gnomad, .SCHEMA_GNOMAD, "gnomad",
                                     factor_levels = list(
                                       consequence = .CONSEQUENCE_LEVELS)))
  }
  if ("alphamissense" %in% names(obj)) {
    errs <- c(errs, .check_df_schema(obj$alphamissense,
                                     .SCHEMA_ALPHAMISSENSE,
                                     "alphamissense",
                                     factor_levels = list(
                                       am_class = .AM_CLASS_LEVELS)))
  }
  if ("revel" %in% names(obj)) {
    errs <- c(errs, .check_df_schema(obj$revel, .SCHEMA_REVEL, "revel"))
  }
  if ("cadd" %in% names(obj)) {
    errs <- c(errs, .check_df_schema(obj$cadd, .SCHEMA_CADD, "cadd",
                                     factor_levels = list(
                                       consequence = .CONSEQUENCE_LEVELS)))
  }

  # 3. Position bounds. Only check if meta is present and well-formed.
  if (length(errs) == 0L && "meta" %in% names(obj) &&
      "protein_length" %in% names(obj$meta)) {
    L <- obj$meta$protein_length
    for (tbl in c("domains", "clinvar", "gnomad",
                  "alphamissense", "revel", "cadd")) {
      df <- obj[[tbl]]
      if (is.data.frame(df) && nrow(df) > 0L) {
        if (tbl == "domains") {
          if (any(df$start < 1L) || any(df$end > L)) {
            errs <- c(errs, sprintf("%s: positions outside [1, %d]",
                                    tbl, L))
          }
        } else if ("pos" %in% names(df)) {
          if (any(df$pos < 1L) || any(df$pos > L)) {
            errs <- c(errs, sprintf("%s: positions outside [1, %d]",
                                    tbl, L))
          }
        }
      }
    }
  }

  if (length(errs) == 0L) {
    if (strict) return(invisible(TRUE))
    list(valid = TRUE, issues = character(0))
  } else {
    if (strict) {
      stop(paste0("validate_gene_data() failed validation:\n  ",
                  paste(errs, collapse = "\n  ")), call. = FALSE)
    }
    list(valid = FALSE, issues = errs)
  }
}

#' Construct an empty per-gene object conforming to the spec
#'
#' Builds a placeholder list with all 7 required tables present
#' but empty (zero rows, correct columns and types). Useful as a
#' template when writing custom annotation files or as a fallback.
#'
#' @param gene HGNC symbol for the `meta$gene` field.
#' @param uniprot_id UniProt accession.
#' @param protein_length Integer length.
#' @return A list conforming to `DATA_FORMAT_SPEC.md`.
#' @export
empty_gene_data <- function(gene = "UNKNOWN",
                            uniprot_id = NA_character_,
                            protein_length = 0L) {
  list(
    meta = data.frame(
      gene = gene,
      uniprot_id = uniprot_id,
      protein_length = as.integer(protein_length),
      ensembl_gene_id = NA_character_,
      ensembl_transcript_id = NA_character_,
      build_date = Sys.Date(),
      source_versions = "",
      stringsAsFactors = FALSE
    ),
    domains = data.frame(
      start = integer(0), end = integer(0),
      name = character(0), accession = character(0),
      source = factor(character(0), levels = .DOMAIN_SOURCE_LEVELS),
      stringsAsFactors = FALSE
    ),
    clinvar = data.frame(
      pos = integer(0), aa_ref = character(0), aa_alt = character(0),
      aa_change = character(0),
      significance = factor(character(0), levels = .SIGNIFICANCE_LEVELS),
      review_status = factor(character(0), levels = .REVIEW_STATUS_LEVELS),
      clinvar_id = character(0),
      stringsAsFactors = FALSE
    ),
    gnomad = data.frame(
      pos = integer(0), aa_ref = character(0), aa_alt = character(0),
      aa_change = character(0),
      consequence = factor(character(0), levels = .CONSEQUENCE_LEVELS),
      af_joint = numeric(0), ac_joint = integer(0), an_joint = integer(0),
      filter = character(0),
      stringsAsFactors = FALSE
    ),
    alphamissense = data.frame(
      pos = integer(0), aa_ref = character(0), aa_alt = character(0),
      aa_change = character(0),
      am_score = numeric(0),
      am_class = factor(character(0), levels = .AM_CLASS_LEVELS),
      stringsAsFactors = FALSE
    ),
    revel = data.frame(
      pos = integer(0), aa_ref = character(0), aa_alt = character(0),
      aa_change = character(0), revel_score = numeric(0),
      stringsAsFactors = FALSE
    ),
    cadd = data.frame(
      pos = integer(0), aa_ref = character(0), aa_alt = character(0),
      aa_change = character(0),
      consequence = factor(character(0), levels = .CONSEQUENCE_LEVELS),
      cadd_raw = numeric(0), cadd_phred = numeric(0),
      stringsAsFactors = FALSE
    )
  )
}
