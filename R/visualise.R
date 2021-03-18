# Visualise new DMS datasets
# TODO - test

#' Default deepscanscape GGPlot2 theme
#'
#' A GGPlot2 theme based on \code{\link[ggpubr]{theme_pubclean}}, with a generally minimalist style. This theme
#' changes the legend to be on the right, has a central title and subtitle, and removes facet strip and legend
#' key backgrounds.
#'
#' @importFrom ggplot2 %+replace%
#' @export
theme_deepscanscape <- function() {
  ggpubr::theme_pubclean() %+replace%
    ggplot2::theme(legend.position = "right",
                   plot.title = ggplot2::element_text(hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5),
                   strip.background = ggplot2::element_blank(),
                   legend.key = ggplot2::element_blank())
}

#' Plot a heat map of ER scores from a deep mutational scan
#'
#' Plot a heat map showing the ER fitness score for each mutation at each position and the positions mean ER score.
#' The wild type amino acid at each position is marked by outlining the corresponding heatmap tile.
#' This is also the default plot method for \code{\link{deep_mutational_scan}} objects.
#'
#' @param x \code{\link{deep_mutational_scan}}.
#' @return A \code{\link[ggplot2]{ggplot2}} plot.
#' @examples
#' dms <- deepscanscape::deep_scan
#' plot_er_heatmap(dms)
#'
#' # Or equivalently
#' plot(dms)
#'
#' @export
plot_er_heatmap <- function(x) {
  if (!is.deep_mutational_scan(x)) {
    stop("x is not a deep_mutational_scan()")
  }

  df <- tidyr::pivot_longer(x$data[c("position", "wt", amino_acids)],
                            cols = .data$A:.data$Y, names_to = "mut", values_to = "er")
  means <- dplyr::summarise(dplyr::group_by(df, .data$position, , .data$wt),
                            er = mean(.data$er), mut = "Mean", .groups = "drop")
  df <- dplyr::bind_rows(df, means)
  df$mut <- factor(df$mut, levels = c(amino_acids, "", "Mean"))

  limit <- rep(max(abs(df$er)), 2) * c(-1, 1)
  ggplot2::ggplot(mapping = ggplot2::aes(x = .data$position, y = as.integer(.data$mut), fill = .data$er)) +
    ggplot2::geom_tile(data = df[df$wt != df$mut, ]) +
    ggplot2::geom_tile(data = df[df$wt == df$mut, ], mapping = ggplot2::aes(colour = "WT")) +
    ggplot2::scale_fill_distiller(name = "ER", limits = limit, type = "div", palette = "RdBu", direction = 1) +
    ggplot2::scale_colour_manual(name = "", values = c(WT = "black")) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(fill = "white"))) +
    ggplot2::scale_y_continuous(breaks = 1:22, labels = c(amino_acids, "", "Mean")) +
    ggplot2::labs(title = x$gene, x = "Position", y = "") +
    ggplot2::xlim(min(df$position), max(df$position)) +
    theme_deepscanscape() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank())
}

# TODO special density feature for point density?
#' Plot a new study on the deep mutational landscape
#'
#' Project new data onto the deep mutational landscape derived from the combined data, viewed in UMAP space. This view
#' allows you to determine if new data follows the expected distribution and compare it with the average properties for
#' positions in the same regions of the landscape. This can give insights into the likely properties of positions in the
#' new protein.
#'
#' Points represent positions from the new protein(s), positioned based on predicted UMAP coordinates from the
#' combined landscape model. These are plotted on top of the whole combined landscape, which is represented as grey
#' points when no feature is used and hexbins when one is. These bins are coloured based on the mean value of the
#' chosen feature for all points that fall within the hexagon. If the chosen feature is a FoldX Gibbs Free Energy
#' prediction values are initially clamped to -5 <= x <= 5 because larger values lead to difficult to read scales and
#' are generally thought to represent program quirks rather than real biology.
#'
#' The input data is annotated using \code{\link{annotate}} if it is not already, based on the \code{annotated} flag
#' for deep_mutational_scans and the presence of the required \code{umap1/2} columns for data frames. It is better to
#' save the results of \code{annotate} if doing multiple downstream analyses because can be a relatively slow
#' operation.
#'
#' @param x \code{\link{deep_mutational_scan}} or \link[=rbind.deep_mutational_scan]{combined DMS data frame}.
#' @param name Prefer to identify datasets by study or gene
#' @param feature String name of a numeric feature from the \code{\link{deep_mutational_scan}} dataset to project onto
#' the background landscape.
#' @return A \code{\link[ggplot2]{ggplot2}} plot.
#' @examples
#' dms <- deepscanscape::deep_scan
#'
#' # Plot point positions
#' plot_landscape(dms, name = "gene")
#'
#' # Plot against mean conservation
#' dms <- annotate(dms)
#' plot_landscape(dms, feature = "mean_sift")
#'
#' @export
# TODO add examples with multiple studies when more example datasets available
plot_landscape <- function(x, name = c("study", "gene"), feature = NULL) {
  if (is.deep_mutational_scan(x)) {
    if (!x$annotated) {
      warning("deep_mutational_scan is not annotated. Annotating using annotate().")
      x <- annotate(x)
    }
    df <- tibble::as_tibble(x, full = TRUE)
  } else if (is.data.frame(x)) {
    if (!all(c("umap1", "umap2") %in% names(x))) {
      warning("data frame is not annotated. Annotating using annotate().")
      x <- annotate(x)
    }
    df <- validate_combined_dms(x)
  }

  name <- match.arg(name)
  df$name <- df[[name]]
  df$name[is.na(df$name)] <- df[[grep(name, c("study", "gene"), value = TRUE, invert = TRUE)]][is.na(df$name)]

  n_studies <- length(unique(df$name))
  if (n_studies == 1) {
    col_scale <- ggplot2::scale_colour_manual(name = "", values = "black")
  } else if (n_studies <= 8) {
    col_scale <- ggplot2::scale_colour_brewer(name = "", type = "qual", palette = "Set1")
  } else {
    col_scale <- ggplot2::scale_colour_brewer(name = "", type = "qual", palette = "Paired")
  }

  if (is.null(feature)) {
    p <- ggplot2::ggplot(mapping = ggplot2::aes(x = .data$umap1, y = .data$umap2)) +
      ggplot2::geom_point(data = deepscanscape::deep_mutational_scans, mapping = ggplot2::aes(fill = "Background"),
                          shape = 21, colour = "grey") +
      ggplot2::geom_point(data = df, mapping = ggplot2::aes(colour = .data$name)) +
      col_scale +
      ggplot2::scale_fill_manual(name = "", values = c(Background = "grey")) +
      ggplot2::labs(x = "UMAP1", y = "UMAP2") +
      theme_deepscanscape()
  } else {
    if (!feature %in% names(deepscanscape::deep_mutational_scans)) {
      stop("Feature must be a variable in the deep_mutational_scans dataset")
    }

    if (typeof(deepscanscape::deep_mutational_scans[[feature]]) == "character") {
      stop("Features must be numeric variables (e.g. 'mean_sift' or 'entropy_sidechain')")
    }

    comb_df <- deepscanscape::deep_mutational_scans[!is.na(deepscanscape::deep_mutational_scans[[feature]]), ]
    fill_scale <- get_feature_scale(feature)

    if (feature %in% foldx_terms) {
      comb_df[feature] <- clamp(comb_df[[feature]], -5, 5)
    }

    p <- ggplot2::ggplot(mapping = ggplot2::aes(x = .data$umap1, y = .data$umap2)) +
      ggplot2::stat_summary_hex(data = comb_df, fun = mean, bins = 40, mapping = ggplot2::aes(z = .data[[feature]])) +
      ggplot2::geom_point(data = df, shape = 19, colour = "black") +
      ggplot2::geom_point(data = df, mapping = ggplot2::aes(colour = .data$name), shape = 16) +
      fill_scale +
      col_scale +
      ggplot2::labs(x = "UMAP1", y = "UMAP2") +
      theme_deepscanscape()
  }
  return(p)
}

#' Select colour scales for the features in the combined dataset
#'
#' Select from preconfigured scales for the various features in the combined dataset
#'
#' @param feature Feature to fetch the scale for.
#' @param type Return a fill or a colour scale
#' @keywords internal
get_feature_scale <- function(feature, type = c("fill", "colour")) {
  type <- match.arg(type)
  if (type == "fill") {
    scale_distiller <- ggplot2::scale_fill_distiller
    scale_gradient2 <- ggplot2::scale_fill_gradient2
    scale_gradientn <- ggplot2::scale_fill_gradientn
    scale_continuous <- ggplot2::scale_fill_continuous
  } else if (type == "colour") {
    scale_distiller <- ggplot2::scale_colour_distiller
    scale_gradient2 <- ggplot2::scale_colour_gradient2
    scale_gradientn <- ggplot2::scale_colour_gradientn
    scale_continuous <- ggplot2::scale_colour_continuous
  } else {
    # Should never get here, should be caught be match.arg
    stop("Unrecognised 'type' argument")
  }

  pretty <- stringr::str_to_title(stringr::str_replace_all(feature, "_", " "))

  if (feature == "mean_sift") {
    fill_scale <- scale_distiller(name = expression(log[10] * "SIFT"), type = "seq", palette = "RdPu")
  } else if (feature == "mean_score") {
    fill_scale <- scale_gradient2(name = "Mean ER", low = "#d73027", mid = "#ffffbf", high = "#4575b4")
  } else if (feature %in% amino_acids) {
    fill_scale <- scale_gradient2(name = feature, low = "#4575b4", mid = "#ffffbf", high = "#d73027")
  } else if (feature == "all_atom_rel") {
    fill_scale <- scale_gradientn(colours = c("#1a2a6c", "#b21f1f", "#fdbb2d"), values = c(0, 0.4, 1))
  } else if (feature %in% foldx_terms) {
    fill_scale <- scale_gradient2(name = pretty, low = "#4575b4", mid = "#ffffbf", high = "#d73027")
  } else if (feature == "hydrophobicity") {
    fill_scale <- scale_gradientn(colours = c("#4575b4", "#e0f3f8", "#fee090", "#fc8d59", "#d73027"),
                                  values = scales::rescale(c(-0.4, 0, 0.4, 0.8, 1.2)),
                                  limits = c(-0.4, 1.201))
  } else if (feature %in% c("phi",  "psi")) {
    fill_scale <- scale_gradient2(name = pretty, low = "#c51b7d", mid = "#ffffbf", high = "#4d9221")
  } else {
    fill_scale <- scale_continuous(name = pretty)
  }

  return(fill_scale)
}

#' Plot amino acid subtype frequencies
#'
#' Plot the frequency each cluster occurs in an annotated deep mutational scan dataset, in order to identify data that
#' strays from the expected frequencies. For example large proportions of outliers could suggest abnormal data or many
#' permissive positions a weakly conserved protein.
#'
#' @param x \code{\link{deep_mutational_scan}} or \link[=rbind.deep_mutational_scan]{combined DMS data frame}.
#' @return A \code{\link[ggplot2]{ggplot2}} plot.
#' @examples
#' dms <- annotate(deepscanscape::deep_scan)
#' plot_cluster_frequencies(dms)
#'
#' @export
# TODO - add examples using multiple datasets when examples available
plot_cluster_frequencies <- function(x) {
  if (is.deep_mutational_scan(x)) {
    if (!x$annotated) {
      warning("deep_mutational_scan is not annotated. Annotating using annotate().")
      x <- annotate(x)
    }
    df <- x$data
  } else if (is.data.frame(x)) {
    df <- validate_combined_dms(x, annotated = TRUE)
  }

  overview <- dplyr::summarise(dplyr::group_by(df, .data$wt, .data$cluster), n = dplyr::n(), .groups = "drop_last")
  overview <- dplyr::ungroup(dplyr::mutate(overview, prop = .data$n / sum(.data$n)))
  overview$wt <- factor(overview$wt, levels = sort(unique(overview$wt), decreasing = TRUE))
  overview$cluster_num <- stringr::str_sub(overview$cluster, start = 2)

  max_cluster <- max(as.integer(stringr::str_subset(overview$cluster_num, "[0-9]+")))
  overview$cluster_num <- factor(overview$cluster_num, levels = c("O", "A", "P", rev(seq_len(max_cluster))))

  counts <- dplyr::summarise(dplyr::group_by(overview, .data$wt), n = sum(.data$n), .groups = "drop")
  sec_axis <- ggplot2::dup_axis(name = "", labels = counts$n)

  if (max_cluster <= 8) {
    fill_scale <- ggplot2::scale_fill_manual(values = c("1" = "#e41a1c", "2" = "#377eb8", "3" = "#4daf4a",
                                                        "4" = "#984ea3", "5" = "#ff7f00", "6" = "#ffff33",
                                                        "7" = "#42b7ce", "8" = "#f781bf",
                                                        "P" = "#adadad", "A" = "#a65628", "O" = "#666666"))
  } else {
    fill_scale <- ggplot2::scale_fill_brewer(name = "Subtype")
  }

  ggplot2::ggplot(overview, ggplot2::aes(x = as.integer(.data$wt), y = .data$prop, fill = .data$cluster_num)) +
    ggplot2::scale_x_continuous(breaks = seq_len(20), labels = counts$wt, sec.axis = sec_axis) +
    ggplot2::scale_y_continuous(labels = function(x) stringr::str_c(100 * x, "%")) +
    ggplot2::coord_flip(expand = FALSE) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::labs(x = "WT Amino Acid", y = "Percentage of Positions") +
    ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE)) +
    fill_scale +
    theme_deepscanscape() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank())
}
