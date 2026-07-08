## ===================================================================
## msaVariant: plot_variant_overlay()
## ===================================================================
##
## Direct port of the verified make_publication_figure.R into a
## package function. The body of this function is the script verbatim;
## only the gene-specific inputs (gene symbol, aligned FASTA, ref
## sequence name, variant position+label, window) and the domain
## display logic are parameterized.
##
## DESIGN PRINCIPLE (from the working script):
##   1. Every panel—including the MSA—uses standard ggplot2 primitives
##      sharing the exact same continuous x-axis coordinate system.
##   2. All text labels are isolated into a dedicated left-hand column
##      so patchwork's matrix layout aligns labels to data perfectly.
##   3. Legend dashboard is constructed inside the same matrix so it
##      cannot drift outside the data columns.

#' Generate a publication-grade variant-overlay figure
#'
#' Composes a multi-panel figure (variant header, MSA, annotation
#' tracks, legend dashboard) for a clinical variant in the context
#' of evolutionary conservation and population/clinical annotations.
#' The gene's annotation data is read from the package cache via
#' [fetch_gene_data()]. The user supplies only the gene symbol, the
#' aligned FASTA, and the variant position.
#'
#' @param gene Gene symbol (string). Must match the name of a `.rds`
#'   bundle in the cache. Use [import_local_bundle()] to add a local
#'   bundle, or [fetch_gene_data()] to download from Zenodo.
#' @param aligned_fasta Path to an aligned FASTA file (e.g. MUSCLE
#'   output). The reference sequence (matching the canonical UniProt
#'   isoform) should be one of the sequences in the file.
#' @param variant_pos Residue position in the reference (UniProt
#'   canonical numbering). Can also be a string like `"R175H"` or
#'   `"p.R175H"` from which the position is parsed.
#' @param variant_label Optional display label (e.g. `"p.R175H"`).
#'   If `NULL`, derived from `variant_pos` and the reference residue.
#' @param aa_change Optional variant identifier used for the ACMG
#'   evidence strip, e.g. `"R175H"`. If `NULL`, it is parsed from
#'   `variant_label` when that carries a full ref+pos+alt substitution.
#'   ACMG codes that require the alternate allele are only shown when a
#'   complete `aa_change` is available.
#' @param pp3_min_predictors Integer 1–3 (default 2) forwarded to
#'   [compute_acmg_codes()]: PP3 fires when at least this many of the
#'   three PP3 predictors (AlphaMissense/REVEL/CADD) pass their
#'   threshold.
#' @param color_scheme Figure color preset: `"journal"` (default,
#'   publication-quality, colorblind-safe, grayscale-legible, muted),
#'   `"colorblind"` (explicit Okabe-Ito), or `"grayscale"` (pure
#'   black-and-white for print).
#' @param acmg_colors Optional named vector overriding ACMG chip colors
#'   (names among `PS1`,`PM1`,`PM2`,`PM5`,`PP3`); partial allowed.
#'   `NULL` uses the scheme.
#' @param variant_highlight_color Optional single color for the variant
#'   highlight fill; `NULL` uses the scheme.
#' @param track_palettes Optional named list overriding individual track
#'   palettes: any of `clinvar`, `gnomad`, `cadd`, `alphamissense`,
#'   `revel`. `NULL` uses the scheme.
#' @param msa_color_scheme Optional MSA coloring: a scheme name
#'   (`"journal"`, `"colorblind"`, `"clustal"`, `"grayscale"`) or an
#'   explicit named per-residue color vector. `NULL` uses the scheme.
#' @param show_acmg,show_domain,show_clinvar,show_gnomad,show_cadd,show_alphamissense,show_revel
#'   Logical toggles (all `TRUE` by default) for each annotation layer.
#'   Setting one to `FALSE` fully removes that layer and the figure
#'   re-flows (patchwork heights and row assignments adjust to the
#'   number of active layers — no blank gaps). The variant header and
#'   the MSA are always shown.
#' @param acmg_codes Character vector of ACMG codes to *display*
#'   (default all five: `c("PS1","PM1","PM2","PM5","PP3")`). The full
#'   ACMG computation always runs internally; only codes in this set
#'   that were actually triggered are drawn as chips.
#' @param window Length-2 integer vector `c(start, end)` defining
#'   the residue range to display. Defaults to ±20 residues around
#'   `variant_pos`, clipped to the protein.
#' @param ref_name Name of the reference sequence in the FASTA.
#'   If `NULL`, the sequence whose name contains the gene symbol is
#'   used; otherwise the first sequence.
#' @param title Plot title.
#' @param subtitle Plot subtitle.
#'
#' @return A `patchwork` plot object. Save with `ggplot2::ggsave()`.
#'
#' @examples
#' ## Runnable with the synthetic DEMO1 example bundle shipped in the
#' ## package (fabricated data; not real predictions). Use a temporary
#' ## cache so the example does not touch your real cache directory.
#' Sys.setenv(MSAVARIANT_CACHE = tempfile("msaVariant_cache_"))
#' import_local_bundle(
#'   system.file("extdata", "DEMO1.rds", package = "msaVariant"),
#'   gene = "DEMO1"
#' )
#' p <- plot_variant_overlay(
#'   gene          = "DEMO1",
#'   aligned_fasta = system.file("extdata", "demo_aligned.fasta",
#'                               package = "msaVariant"),
#'   variant_pos   = 21,
#'   variant_label = "p.R21H"
#' )
#'
#' \dontrun{
#' ## Save a publication-resolution figure
#' ggplot2::ggsave("demo1_R21H.png", p,
#'                 width = 15, height = 9, dpi = 300, bg = "white")
#' }
#'
#' @export
plot_variant_overlay <- function(gene,
                                  aligned_fasta,
                                  variant_pos,
                                  variant_label = NULL,
                                  window        = NULL,
                                  ref_name      = NULL,
                                  title         = NULL,
                                  subtitle      = NULL,
                                  aa_change     = NULL,
                                  pp3_min_predictors = 2L,
                                  color_scheme  = "journal",
                                  acmg_colors   = NULL,
                                  variant_highlight_color = NULL,
                                  track_palettes = NULL,
                                  msa_color_scheme = NULL,
                                  show_acmg          = TRUE,
                                  show_domain        = TRUE,
                                  show_clinvar       = TRUE,
                                  show_gnomad        = TRUE,
                                  show_cadd          = TRUE,
                                  show_alphamissense = TRUE,
                                  show_revel         = TRUE,
                                  acmg_codes = c("PS1","PM1","PM2","PM5","PP3")) {
  ## --- Input parsing & defaults ---
  stopifnot(is.character(gene), length(gene) == 1L)
  stopifnot(file.exists(aligned_fasta))

  ## Resolve the color specification (preset + optional overrides).
  ## Presentation only — never affects ACMG logic.
  col <- resolve_plot_colors(
    color_scheme            = color_scheme,
    acmg_colors             = acmg_colors,
    variant_highlight_color = variant_highlight_color,
    track_palettes          = track_palettes,
    msa_color_scheme        = msa_color_scheme
  )

  ## Parse variant_pos if given as a string
  if (!is.numeric(variant_pos)) {
    s <- sub("^p\\.", "", as.character(variant_pos))
    m <- regmatches(s, regexec("([0-9]+)", s))[[1]]
    if (length(m) < 1) stop("Could not parse variant position from: ", variant_pos)
    variant_pos <- as.integer(m[1])
  }
  variant_pos <- as.integer(variant_pos)
  
  ## Load the gene bundle from cache
  bundle <- fetch_gene_data(gene)
  protein_length <- bundle$meta$protein_length[1]
  if (variant_pos < 1L || variant_pos > protein_length) {
    stop(sprintf("variant_pos %d outside protein (length %d)",
                 variant_pos, protein_length))
  }
  
  ## Resolve reference sequence in the FASTA
  aa <- Biostrings::readAAStringSet(aligned_fasta)
  seq_names <- names(aa)
  if (is.null(ref_name)) {
    hits <- grep(toupper(gene), toupper(seq_names), fixed = TRUE)
    ref_name <- if (length(hits)) seq_names[hits[1]] else seq_names[1]
  }
  if (!(ref_name %in% seq_names)) {
    stop(sprintf("ref_name '%s' not found in FASTA. Sequences: %s",
                 ref_name, paste(seq_names, collapse = ", ")))
  }
  
  ## Window defaults to ±20 around variant
  if (is.null(window)) {
    window <- c(max(1L, variant_pos - 20L),
                min(as.integer(protein_length), variant_pos + 20L))
  }
  stopifnot(length(window) == 2L)
  WINDOW_START_RES <- as.integer(window[1])
  WINDOW_END_RES   <- as.integer(window[2])
  VARIANT_POS      <- variant_pos
  
  ## ----- (Script body begins. Verbatim apart from parameterized refs.) -----
  coord_map <- build_msa_coord_map(aligned_fasta, ref_name = ref_name)
  start_col <- min(coord_map$msa_col[coord_map$residue_pos %in% WINDOW_START_RES:WINDOW_END_RES], na.rm = TRUE)
  end_col   <- max(coord_map$msa_col[coord_map$residue_pos %in% WINDOW_START_RES:WINDOW_END_RES], na.rm = TRUE)
  variant_col <- coord_map$msa_col[match(VARIANT_POS, coord_map$residue_pos)]
  
  ## Derive default variant_label if missing
  if (is.null(variant_label)) {
    ref_aa <- substr(as.character(aa[[ref_name]]), variant_col, variant_col)
    VARIANT_LABEL <- sprintf("p.%s%d", ref_aa, VARIANT_POS)
  } else {
    VARIANT_LABEL <- variant_label
  }
  
  ## Read & transform alignment to long form
  n_seqs <- length(aa)
  alignment_list <- lapply(seq_len(n_seqs), function(i) {
    data.frame(
      seq_index = i,
      seq_name  = seq_names[i],
      msa_col   = seq_len(Biostrings::width(aa)[i]),
      residue   = strsplit(as.character(aa[[i]]), "")[[1]],
      stringsAsFactors = FALSE
    )
  })
  df_msa <- do.call(rbind, alignment_list)
  
  df_msa_window <- df_msa[df_msa$msa_col >= start_col & df_msa$msa_col <= end_col, , drop = FALSE]
  df_msa_window <- df_msa_window[!duplicated(df_msa_window[, c("msa_col", "seq_name")]), , drop = FALSE]
  df_msa_window$seq_name <- factor(df_msa_window$seq_name, levels = rev(seq_names))
  
  ## Fetch annotation tables from the bundle
  cv         <- bundle$clinvar
  am         <- bundle$alphamissense
  domains_df <- bundle$domains
  gnomad_df  <- bundle$gnomad
  revel_df   <- bundle$revel
  cadd_df    <- bundle$cadd
  
  ## Shared viewport
  shared_viewport <- ggplot2::coord_cartesian(xlim = c(start_col - 0.5, end_col + 0.5), expand = FALSE)
  
  ## Custom MSA panel
  clustal_colors <- col$msa
  
  p_msa_data <- ggplot2::ggplot(df_msa_window, ggplot2::aes(x = msa_col, y = seq_name)) +
    ggplot2::geom_tile(ggplot2::aes(fill = residue), color = "white", linewidth = 0.1) +
    ggplot2::geom_text(ggplot2::aes(label = residue), family = "mono", fontface = "bold", size = 3) +
    ggplot2::annotate("rect", xmin = variant_col - 0.5, xmax = variant_col + 0.5,
             ymin = 0.5, ymax = n_seqs + 0.5,
             alpha = 0.35, fill = col$variant_fill, color = col$variant_outline, linewidth = 1.2) +
    ggplot2::scale_fill_manual(values = clustal_colors, guide = "none") +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = c(0.5, 0.5))) +
    ggplot2::coord_cartesian(xlim = c(start_col - 0.5, end_col + 0.5),
                    expand = FALSE, clip = "off") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(l = 0, r = 5, t = 2, b = 2)
    )
  
  p_msa_labels <- ggplot2::ggplot(data.frame(
      y = factor(seq_names, levels = rev(seq_names)),
      name = seq_names
    )) +
    ggplot2::geom_text(ggplot2::aes(x = 1, y = y, label = name),
              fontface = "bold", family = "mono", hjust = 1, size = 3.5) +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = c(0.5, 0.5))) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(l = 5, r = 6, t = 2, b = 2))
  
  ## Variant header panel
  p_variant_data <- ggplot2::ggplot() +
    ggplot2::annotate("rect",
             xmin = variant_col - 0.5, xmax = variant_col + 0.5,
             ymin = 0.4, ymax = 1.0,
             alpha = 0.35, fill = col$variant_fill,
             color = col$variant_outline, linewidth = 1.2) +
    ggplot2::annotate("text",
             x = variant_col, y = 1.5,
             label = VARIANT_LABEL,
             color = col$accent, fontface = "bold", size = 5.4) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, 2), expand = c(0, 0)) +
    shared_viewport +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(l = 0, r = 5, t = 4, b = 2))
  
  p_variant_label <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 1, y = 1,
             label = "Variant:",
             fontface = "bold", family = "mono",
             hjust = 1, size = 3.8,
             color = col$accent) +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, 2), expand = c(0, 0)) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(l = 5, r = 6, t = 4, b = 2))
  
  ## ----- ACMG evidence-code strip -----------------------------------
  ## A compact strip directly below the variant header showing the
  ## triggered ACMG codes as colored chips. Codes are computed purely
  ## from the bundle tables via compute_acmg_codes(); untriggered codes
  ## are hidden. This is additive and does not touch any track.

  ## Resolve an aa_change for ACMG evaluation: explicit arg wins, else
  ## parse from the (possibly user-supplied) variant label.
  acmg_aa_change <- if (!is.null(aa_change)) {
    sub("^p\\.", "", as.character(aa_change))
  } else {
    sub("^p\\.", "", VARIANT_LABEL)
  }
  ## Full ACMG computation always runs (logic is independent of display).
  triggered_codes <- tryCatch(
    compute_acmg_codes(bundle, VARIANT_POS, acmg_aa_change,
                       pp3_min_predictors = pp3_min_predictors),
    error = function(e) character(0)
  )
  ## Only the user-requested subset is drawn (preserve canonical order).
  display_codes <- triggered_codes[triggered_codes %in% acmg_codes]

  ## Evidence-strength palette comes from the resolved color scheme
  ## (strong = darkest/warmest, supporting = light amber).
  acmg_palette <- col$acmg

  p_acmg_label <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 1, y = 1,
             label = "ACMG:",
             fontface = "bold", family = "mono",
             hjust = 1, size = 3.8, color = "#37474F") +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, 2), expand = c(0, 0)) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(l = 5, r = 6, t = 2, b = 2))

  if (length(display_codes) > 0L) {
    chip_df <- data.frame(
      code = factor(display_codes, levels = display_codes),
      x    = seq_along(display_codes),
      fill = unname(acmg_palette[display_codes]),
      txt  = unname(.contrast_text(acmg_palette[display_codes])),
      stringsAsFactors = FALSE
    )
    ## Left-align chips within the data column regardless of count.
    x_hi <- max(9, length(display_codes) + 0.5)
    p_acmg_data <- ggplot2::ggplot(chip_df) +
      ggplot2::geom_tile(ggplot2::aes(x = x, y = 1, fill = code),
               width = 0.86, height = 0.7) +
      ggplot2::geom_text(ggplot2::aes(x = x, y = 1, label = code, color = txt),
               fontface = "bold", size = 3.6) +
      ggplot2::scale_color_identity() +
      ggplot2::scale_fill_manual(values = stats::setNames(chip_df$fill, chip_df$code),
               guide = "none") +
      ggplot2::scale_x_continuous(limits = c(0.4, x_hi), expand = c(0, 0)) +
      ggplot2::scale_y_continuous(limits = c(0.5, 1.5), expand = c(0, 0)) +
      ggplot2::theme_void() +
      ggplot2::theme(plot.margin = ggplot2::margin(l = 0, r = 5, t = 2, b = 2))
  } else {
    ## No codes triggered: an empty strip with a muted placeholder so
    ## the row does not read as a rendering failure.
    p_acmg_data <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 1,
               label = "(no ACMG codes triggered)",
               hjust = 0, size = 3.2, color = "grey55", fontface = "italic") +
      ggplot2::scale_x_continuous(limits = c(0.4, 9), expand = c(0, 0)) +
      ggplot2::scale_y_continuous(limits = c(0.5, 1.5), expand = c(0, 0)) +
      ggplot2::theme_void() +
      ggplot2::theme(plot.margin = ggplot2::margin(l = 0, r = 5, t = 2, b = 2))
  }

  ## Track factory (closes over aligned_fasta, ref_name, variant_col, shared_viewport)
  make_track_panel <- function(df, value_col, track_title, custom_scale) {
    if (is.null(df) || nrow(df) == 0L) return(list(label = NULL, data = NULL))
    df$msa_col <- map_variant_to_msa(df$pos, aligned_fasta, ref_name)
    
    p_label <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 1, y = 1, label = track_title, fontface = "bold", hjust = 0, size = 3.8) +
      ggplot2::theme_void() +
      ggplot2::theme(plot.margin = ggplot2::margin(l = 5, r = 10, t = 0, b = 0))
    
    p_data <- ggplot2::ggplot(df) +
      ggplot2::geom_tile(ggplot2::aes(x = msa_col, y = 1, fill = .data[[value_col]]), height = 0.8) +
      ggplot2::annotate("rect", xmin = variant_col - 0.5, xmax = variant_col + 0.5,
               ymin = 0.5, ymax = 1.5, alpha = 0.15, fill = col$variant_outline, color = col$variant_outline, linewidth = 0.3) +
      custom_scale +
      shared_viewport +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(l = 0, r = 5, t = 2, b = 2)
      )
    return(list(label = p_label, data = p_data))
  }
  
  ## Domain panel — pick the most-informative domain that overlaps
  ## the visible window (generalization of the original "use 3rd domain"
  ## hardcoding, which was TP53-specific).
  if (!is.null(domains_df) && nrow(domains_df) > 0L) {
    domains_df$msa_start <- map_variant_to_msa(domains_df$start, aligned_fasta, ref_name)
    domains_df$msa_end   <- map_variant_to_msa(domains_df$end,   aligned_fasta, ref_name)
    
    visible <- !is.na(domains_df$msa_start) & !is.na(domains_df$msa_end) &
               domains_df$msa_end >= start_col & domains_df$msa_start <= end_col
    visible_doms <- domains_df[visible, , drop = FALSE]
    
    if (nrow(visible_doms) > 0L) {
      ## Prefer Pfam, then domains containing the variant, then the
      ## largest (most specific) overlap with the window.
      is_pfam <- !is.na(visible_doms$source) &
                 as.character(visible_doms$source) == "Pfam"
      is_super <- grepl("superfamily|repeat-like|protein\\.[0-9]+|containing.protein",
                         visible_doms$name, ignore.case = TRUE)
      candidate <- visible_doms[(is_pfam | !is_super), , drop = FALSE]
      if (nrow(candidate) == 0L) candidate <- visible_doms
      contains_var <- candidate$msa_start <= variant_col &
                      candidate$msa_end   >= variant_col
      dom_show <- if (any(contains_var)) {
        candidate[contains_var, , drop = FALSE][1, , drop = FALSE]
      } else {
        candidate[1, , drop = FALSE]
      }
      visible_start <- max(dom_show$msa_start, start_col)
      visible_end   <- min(dom_show$msa_end, end_col)
      text_center   <- (visible_start + visible_end) / 2
      DOMAIN_DISPLAY_LABEL <- dom_show$name
    } else {
      visible_doms <- domains_df[1, , drop = FALSE]
      visible_start <- start_col; visible_end <- end_col
      text_center   <- (start_col + end_col) / 2
      DOMAIN_DISPLAY_LABEL <- "(no domain in window)"
    }
    
    p_dom_label <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 1, y = 1, label = "Domains:", fontface = "bold", hjust = 0, size = 3.8) +
      ggplot2::theme_void() +
      ggplot2::theme(plot.margin = ggplot2::margin(l = 5, r = 10, t = 0, b = 0))
    
    p_dom_data <- ggplot2::ggplot(visible_doms) +
      ggplot2::geom_rect(ggplot2::aes(xmin = msa_start - 0.5, xmax = msa_end + 0.5, ymin = 0.6, ymax = 1.4),
                fill = col$domain, color = "black", linewidth = 0.4) +
      ggplot2::annotate("text", x = text_center, y = 1.0, label = DOMAIN_DISPLAY_LABEL,
               color = "black", fontface = "bold.italic", size = 3.5) +
      ggplot2::annotate("rect", xmin = variant_col - 0.5, xmax = variant_col + 0.5,
               ymin = 0.5, ymax = 1.5, alpha = 0.15, fill = col$variant_outline, color = col$variant_outline, linewidth = 0.3) +
      shared_viewport +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(l = 0, r = 5, t = 2, b = 2)
      )
  } else {
    p_dom_label <- NULL; p_dom_data <- NULL
  }
  
  ## ClinVar track
  cv_palette <- col$clinvar
  cv_pair <- make_track_panel(cv, "significance", "ClinVar",
                              ggplot2::scale_fill_manual(values = cv_palette,
                                                          name = "ClinVar"))
  
  ## gnomAD track
  if (!is.null(gnomad_df) && "af_joint" %in% names(gnomad_df)) {
    gnomad_agg <- aggregate(af_joint ~ pos, gnomad_df, max)
    gnomad_agg$af_log <- log10(gnomad_agg$af_joint)
    gnomad_pair <- make_track_panel(gnomad_agg, "af_log", "gnomAD",
                                    ggplot2::scale_fill_gradientn(
                                      colors = col$gnomad,
                                      name = "gnomAD log10(AF)"))
  } else {
    gnomad_pair <- list(label = NULL, data = NULL)
  }
  
  ## CADD track
  if (!is.null(cadd_df)) {
    cadd_agg <- aggregate(cadd_phred ~ pos, cadd_df, median)
    cadd_pair <- make_track_panel(cadd_agg, "cadd_phred", "CADD",
                                  ggplot2::scale_fill_gradientn(
                                    colors = col$cadd,
                                    limits = c(0,40), name = "CADD PHRED"))
  } else {
    cadd_pair <- list(label = NULL, data = NULL)
  }
  
  ## AlphaMissense track
  if (!is.null(am)) {
    am_agg <- aggregate(am_score ~ pos, am, mean)
    am_pair <- make_track_panel(am_agg, "am_score", "AlphaMissense",
                                ggplot2::scale_fill_gradientn(
                                  colors = col$alphamissense,
                                  limits = c(0,1), name = "AlphaMissense"))
  } else {
    am_pair <- list(label = NULL, data = NULL)
  }
  
  ## REVEL track
  if (!is.null(revel_df)) {
    revel_agg <- aggregate(revel_score ~ pos, revel_df, mean)
    revel_pair <- make_track_panel(revel_agg, "revel_score", "REVEL",
                                   ggplot2::scale_fill_gradientn(
                                     colors = col$revel,
                                     limits = c(0,1), name = "REVEL"))
  } else {
    revel_pair <- list(label = NULL, data = NULL)
  }
  
  ## ---------------------------------------------------------------
  ## Assemble only the ACTIVE layers so the figure re-flows cleanly.
  ## A track row is active iff its toggle is TRUE and the bundle has
  ## data for it. The variant header and MSA are always shown; the
  ## ACMG strip follows show_acmg.
  ## ---------------------------------------------------------------
  track_rows <- list()
  add_track <- function(active, label, data, key) {
    if (isTRUE(active) && !is.null(data))
      track_rows[[length(track_rows) + 1L]] <<-
        list(key = key, label = label, data = data)
  }
  add_track(show_domain,        p_dom_label,       p_dom_data,       "domain")
  add_track(show_clinvar,       cv_pair$label,     cv_pair$data,     "clinvar")
  add_track(show_gnomad,        gnomad_pair$label, gnomad_pair$data, "gnomad")
  add_track(show_cadd,          cadd_pair$label,   cadd_pair$data,   "cadd")
  add_track(show_alphamissense, am_pair$label,     am_pair$data,     "alphamissense")
  add_track(show_revel,         revel_pair$label,  revel_pair$data,  "revel")

  ## Build one color legend per active track (domains have none), from
  ## the still-legended data panel, then strip legends off the panels.
  strip_legend <- function(p) if (!is.null(p)) p + ggplot2::theme(legend.position = "none") else NULL
  active_legends <- list()
  for (i in seq_along(track_rows)) {
    tr <- track_rows[[i]]
    if (!identical(tr$key, "domain")) {
      active_legends[[length(active_legends) + 1L]] <-
        cowplot::get_legend(tr$data + ggplot2::theme(legend.direction = "horizontal"))
    }
    track_rows[[i]]$data <- strip_legend(tr$data)
  }

  ## Residue-position x-axis goes on the lowest active track, or on the
  ## MSA panel when no tracks are shown, so it is never left orphaned.
  add_axis <- function(p) {
    pretty_pos <- pretty(c(WINDOW_START_RES, WINDOW_END_RES), n = 5)
    pretty_pos <- pretty_pos[pretty_pos >= WINDOW_START_RES & pretty_pos <= WINDOW_END_RES]
    tick_cols  <- coord_map$msa_col[match(pretty_pos, coord_map$residue_pos)]
    keep <- !is.na(tick_cols); pretty_pos <- pretty_pos[keep]; tick_cols <- tick_cols[keep]
    p +
      ggplot2::scale_x_continuous(breaks = tick_cols, labels = pretty_pos, expand = c(0, 0)) +
      ggplot2::theme(
        axis.text.x  = ggplot2::element_text(size = 9, face = "bold", color = "black"),
        axis.title.x = ggplot2::element_text(size = 11, face = "bold", color = "black",
                                             margin = ggplot2::margin(t = 8)),
        axis.ticks.x = ggplot2::element_line(color = "black")
      ) +
      ggplot2::xlab(sprintf("%s residue position", gene))
  }
  if (length(track_rows) > 0L) {
    li <- length(track_rows)
    track_rows[[li]]$data <- add_axis(track_rows[[li]]$data)
  } else {
    p_msa_data <- add_axis(p_msa_data)
  }

  ## Legend dashboard — only for the active tracks that carry a legend.
  have_legends <- length(active_legends) > 0L
  if (have_legends) {
    p_leg_row1 <- cowplot::plot_grid(
      plotlist = active_legends[seq_len(min(3, length(active_legends)))],
      nrow = 1, align = "h", rel_widths = rep(1, min(3, length(active_legends))))
    if (length(active_legends) > 3) {
      p_leg_row2 <- cowplot::plot_grid(
        plotlist = active_legends[4:length(active_legends)],
        nrow = 1, align = "h", rel_widths = rep(1, length(active_legends) - 3))
      p_panel3_data <- cowplot::plot_grid(p_leg_row1, p_leg_row2, ncol = 1, rel_heights = c(1, 1))
    } else {
      p_panel3_data <- p_leg_row1
    }
    p_panel3_label <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 1, y = 1, label = "Color Legends:", fontface = "bold", hjust = 0, size = 3.8) +
      ggplot2::theme_void() +
      ggplot2::theme(plot.margin = ggplot2::margin(l = 5, r = 10, t = 0, b = 0))
  }

  ## Ordered rows: variant (always), ACMG (optional), MSA (always),
  ## active tracks, legends (if any). Heights + row letters follow the
  ## actual number of active rows, so nothing leaves a blank gap.
  rows <- list()
  add_row <- function(label, data, height)
    rows[[length(rows) + 1L]] <<- list(label = label, data = data, height = height)

  add_row(p_variant_label, p_variant_data, 0.7)
  if (isTRUE(show_acmg)) add_row(p_acmg_label, p_acmg_data, 0.55)
  add_row(p_msa_labels, p_msa_data, 4.5)
  for (tr in track_rows) add_row(tr$label, tr$data, 0.5)
  if (have_legends) add_row(p_panel3_label, p_panel3_data, 1.8)

  active_components <- unlist(lapply(rows, function(r) list(r$label, r$data)),
                              recursive = FALSE)
  row_heights <- vapply(rows, function(r) r$height, numeric(1))
  n_active_rows <- length(rows)

  layout_chars <- LETTERS
  design_string <- ""
  char_idx <- 1
  for (r in seq_len(n_active_rows)) {
    design_string <- paste0(design_string, layout_chars[char_idx], layout_chars[char_idx + 1], "\n")
    char_idx <- char_idx + 2
  }
  
  ## Title/subtitle defaults
  if (is.null(title)) {
    title <- sprintf("%s %s: Clinical-Genetics Annotation Overlay Stack",
                     gene, VARIANT_LABEL)
  }
  if (is.null(subtitle)) {
    subtitle <- sprintf("%d-species MSA window, residues %d..%d",
                        n_seqs, WINDOW_START_RES, WINDOW_END_RES)
  }
  
  final_publication_plot <- patchwork::wrap_plots(
      active_components,
      design  = design_string,
      heights = row_heights,
      widths  = c(1.8, 10)
    ) +
    patchwork::plot_annotation(
      title    = title,
      subtitle = subtitle,
      theme = ggplot2::theme(
        plot.title    = ggplot2::element_text(face = "bold", size = 14,
                                              margin = ggplot2::margin(b = 4, l = 5)),
        plot.subtitle = ggplot2::element_text(size = 11, color = "grey30",
                                              margin = ggplot2::margin(b = 10, l = 5))
      )
    )
  
  final_publication_plot
}
