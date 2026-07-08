## ===================================================================
## msaVariant: publication-quality color system
## ===================================================================
##
## resolve_plot_colors() returns a fully-specified color list used by
## plot_variant_overlay(). Three presets are provided; per-element
## overrides take precedence over the preset. Defaults are chosen to be:
##   * colorblind-safe (Okabe-Ito informed),
##   * legible in grayscale (ACMG chips & sequential ramps are
##     luminance-monotonic, so they survive black-and-white printing),
##   * muted / professional (no neon), and
##   * semantically intuitive for ACMG (stronger evidence = darker,
##     warmer red; moderate = orange; supporting = amber).
##
## Nothing here affects ACMG logic — color is presentation only.

## Okabe-Ito reference palette (for documentation / colorblind preset).
.OKABE_ITO <- c(
  black   = "#000000", orange = "#E69F00", skyblue = "#56B4E9",
  green   = "#009E73", yellow = "#F0E442", blue    = "#0072B2",
  vermillion = "#D55E00", purple = "#CC79A7"
)

## Amino-acid property groups (shared by all MSA palettes).
.AA_GROUPS <- list(
  hydrophobic = c("A","I","L","M","F","W","V"),
  positive    = c("K","R"),
  negative    = c("D","E"),
  polar       = c("N","Q","S","T"),
  cysteine    = c("C"),
  glycine     = c("G"),
  proline     = c("P"),
  aromatic    = c("H","Y")
)

## Expand a per-group color vector to a per-residue named vector.
.expand_aa <- function(group_cols, gap = "#FFFFFF") {
  out <- character(0)
  for (g in names(.AA_GROUPS)) {
    for (aa in .AA_GROUPS[[g]]) out[aa] <- group_cols[[g]]
  }
  out["-"] <- gap
  out
}

## Built-in MSA palettes keyed by scheme/name.
.msa_palette <- function(name) {
  name <- tolower(name)
  if (name %in% c("journal", "muted")) {
    return(.expand_aa(list(
      hydrophobic = "#8AA9C9", positive = "#C97B7B", negative = "#B98FB9",
      polar = "#8FBC8F", cysteine = "#D8C57A", glycine = "#E0A96D",
      proline = "#D6D67A", aromatic = "#7FB0B0")))
  }
  if (name %in% c("colorblind", "okabe", "okabe-ito")) {
    return(.expand_aa(list(
      hydrophobic = "#56B4E9", positive = "#D55E00", negative = "#CC79A7",
      polar = "#009E73", cysteine = "#F0E442", glycine = "#E69F00",
      proline = "#7F7F7F", aromatic = "#0072B2")))
  }
  if (name %in% c("clustal", "classic")) {
    return(c(
      "A"="#80a0f0","R"="#f01505","N"="#00a500","D"="#c048c0","C"="#ffff00",
      "Q"="#00a500","E"="#c048c0","G"="#f09048","H"="#15a4a4","I"="#80a0f0",
      "L"="#80a0f0","K"="#f01505","M"="#80a0f0","F"="#80a0f0","P"="#ffff00",
      "S"="#00a500","T"="#00a500","W"="#80a0f0","Y"="#15a4a4","V"="#80a0f0",
      "-"="#ffffff"))
  }
  if (name %in% c("grayscale", "greyscale", "grey", "gray", "bw")) {
    ## Pure B&W: rely on the bold letter, not color. Uniform light tile.
    return(.expand_aa(list(
      hydrophobic = "#F0F0F0", positive = "#F0F0F0", negative = "#F0F0F0",
      polar = "#F0F0F0", cysteine = "#F0F0F0", glycine = "#F0F0F0",
      proline = "#F0F0F0", aromatic = "#F0F0F0"), gap = "#FFFFFF"))
  }
  stop("Unknown msa_color_scheme: ", name)
}

## Full preset definitions.
.scheme_defs <- function(scheme) {
  scheme <- tolower(scheme)
  if (scheme == "journal") {
    list(
      acmg = c(PS1 = "#67000D", PM1 = "#A50026", PM2 = "#D73027",
               PM5 = "#F46D43", PP3 = "#FDAE61"),
      variant_fill = "#E69F00", variant_outline = "#7F0000",
      accent = "#7F0000", domain = "#8FBBA9",
      msa = "journal",
      clinvar = c(Pathogenic="#B2182B", Likely_pathogenic="#EF8A62",
                  VUS="#9E9E9E", Likely_benign="#67A9CF", Benign="#2166AC",
                  Conflicting="#C7A76C", Other="#D9D9D9"),
      gnomad = c("#FCFBFD","#DADAEB","#9E9AC8","#6A51A3","#3F007D"),
      cadd   = c("#F7FCF5","#C7E9C0","#74C476","#31A354","#006D2C"),
      alphamissense = c("#FFF5F0","#FCBBA1","#FB6A4A","#CB181D","#67000D"),
      revel         = c("#FFF5F0","#FCBBA1","#FB6A4A","#CB181D","#67000D")
    )
  } else if (scheme == "colorblind") {
    list(
      acmg = c(PS1 = "#000000", PM1 = "#D55E00", PM2 = "#CC79A7",
               PM5 = "#E69F00", PP3 = "#F0E442"),
      variant_fill = "#E69F00", variant_outline = "#D55E00",
      accent = "#D55E00", domain = "#009E73",
      msa = "colorblind",
      clinvar = c(Pathogenic="#D55E00", Likely_pathogenic="#E69F00",
                  VUS="#999999", Likely_benign="#56B4E9", Benign="#0072B2",
                  Conflicting="#CC79A7", Other="#BBBBBB"),
      gnomad = c("#F7FBFF","#C6DBEF","#6BAED6","#2171B5","#08306B"),
      cadd   = c("#F7FCF5","#C7E9C0","#74C476","#238B45","#00441B"),
      alphamissense = c("#FFF7EC","#FDD49E","#FC8D59","#D7301F","#7F0000"),
      revel         = c("#FFF7EC","#FDD49E","#FC8D59","#D7301F","#7F0000")
    )
  } else if (scheme %in% c("grayscale","greyscale","gray","grey","bw")) {
    list(
      acmg = c(PS1 = "#000000", PM1 = "#3B3B3B", PM2 = "#636363",
               PM5 = "#969696", PP3 = "#CFCFCF"),
      variant_fill = "#BDBDBD", variant_outline = "#000000",
      accent = "#000000", domain = "#D9D9D9",
      msa = "grayscale",
      clinvar = c(Pathogenic="#000000", Likely_pathogenic="#525252",
                  VUS="#969696", Likely_benign="#CCCCCC", Benign="#E8E8E8",
                  Conflicting="#737373", Other="#E0E0E0"),
      gnomad = c("#F0F0F0","#BDBDBD","#737373","#252525","#000000"),
      cadd   = c("#F7F7F7","#CCCCCC","#969696","#525252","#000000"),
      alphamissense = c("#FFFFFF","#CCCCCC","#969696","#525252","#000000"),
      revel         = c("#FFFFFF","#CCCCCC","#969696","#525252","#000000")
    )
  } else {
    stop("Unknown color_scheme: '", scheme,
         "'. Use 'journal', 'colorblind', or 'grayscale'.")
  }
}

#' Resolve the full color specification for a variant-overlay figure
#'
#' @param color_scheme One of `"journal"` (default), `"colorblind"`,
#'   `"grayscale"`.
#' @param acmg_colors Optional named character vector overriding ACMG
#'   chip colors (names among PS1/PM1/PM2/PM5/PP3). Partial is allowed.
#' @param variant_highlight_color Optional single color for the variant
#'   highlight fill.
#' @param track_palettes Optional named list with any of `clinvar`,
#'   `gnomad`, `cadd`, `alphamissense`, `revel`.
#' @param msa_color_scheme Optional MSA palette: a scheme name
#'   (`"journal"`, `"colorblind"`, `"clustal"`, `"grayscale"`) or a
#'   named per-residue color vector.
#' @return A named list of resolved colors.
#' @keywords internal
resolve_plot_colors <- function(color_scheme = "journal",
                                acmg_colors = NULL,
                                variant_highlight_color = NULL,
                                track_palettes = NULL,
                                msa_color_scheme = NULL) {
  col <- .scheme_defs(color_scheme)

  ## MSA palette: preset name, override name, or explicit vector.
  msa_spec <- if (!is.null(msa_color_scheme)) msa_color_scheme else col$msa
  col$msa <- if (is.character(msa_spec) && length(msa_spec) == 1L) {
    .msa_palette(msa_spec)
  } else {
    msa_spec  # explicit named vector
  }

  ## ACMG overrides (merge over preset so partial overrides work).
  if (!is.null(acmg_colors)) {
    for (nm in names(acmg_colors)) col$acmg[[nm]] <- acmg_colors[[nm]]
  }

  ## Variant highlight override (fill only; outline kept from scheme).
  if (!is.null(variant_highlight_color)) col$variant_fill <- variant_highlight_color

  ## Track palette overrides.
  if (!is.null(track_palettes)) {
    for (nm in intersect(names(track_palettes),
                         c("clinvar","gnomad","cadd","alphamissense","revel"))) {
      col[[nm]] <- track_palettes[[nm]]
    }
  }
  col
}

## Choose black or white text for legibility on a given fill color.
## Uses WCAG relative luminance; threshold ~0.5 works well for chips.
.contrast_text <- function(fill) {
  vapply(fill, function(hex) {
    rgb <- grDevices::col2rgb(hex) / 255
    lin <- ifelse(rgb <= 0.03928, rgb/12.92, ((rgb + 0.055)/1.055)^2.4)
    L <- 0.2126*lin[1] + 0.7152*lin[2] + 0.0722*lin[3]
    if (L > 0.5) "#000000" else "#FFFFFF"
  }, character(1))
}
