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
#' @param x A single study \code{\link{deep_mutational_scan}}.
#' @return A \code{\link[ggplot2]{ggplot2}} plot.
#' @examples
#' dms <- deepscanscape::deep_scans$p53
#' plot_er_heatmap(dms)
#'
#' # Or equivalently
#' plot(dms)
#'
#' # For multi_study scans
#' comb_dms <- bind_scans(dms, annotate_missing = TRUE)
#' plot_er_heatmap(dms)
#'
#' @export
plot_er_heatmap <- function(x) {
  if (!is.deep_mutational_scan(x)) {
    stop("x is not a deep_mutational_scan()")
  }

  df <- tidyr::pivot_longer(x$data[c("name", "position", "wt", amino_acids)],
                            cols = .data$A:.data$Y, names_to = "mut", values_to = "er")
  means <- dplyr::summarise(dplyr::group_by(df, .data$name, .data$position, .data$wt),
                            er = mean(.data$er), mut = "Mean", .groups = "drop")
  df <- dplyr::bind_rows(df, means)
  df$mut <- factor(df$mut, levels = c(amino_acids, "", "Mean"))

  limit <- rep(max(abs(df$er)), 2) * c(-1, 1)
  p <- ggplot2::ggplot(mapping = ggplot2::aes(x = .data$position, y = as.integer(.data$mut), fill = .data$er)) +
    ggplot2::geom_tile(data = df[df$wt != df$mut, ]) +
    ggplot2::geom_tile(data = df[df$wt == df$mut, ], mapping = ggplot2::aes(colour = "WT")) +
    ggplot2::scale_fill_distiller(name = "ER", limits = limit, type = "div", palette = "RdBu", direction = 1) +
    ggplot2::scale_colour_manual(name = "", values = c(WT = "black")) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(fill = "white"))) +
    ggplot2::scale_y_continuous(breaks = 1:22, labels = c(amino_acids, "", "Mean")) +
    ggplot2::labs(x = "Position", y = "") +
    theme_deepscanscape() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank())

  if (x$multi_study) {
    p <- p + ggplot2::facet_wrap(~name, scales = "free", ncol = 1)
  } else {
    p <- p + ggplot2::labs(title = x$gene)
  }
  return(p)
}

#' Plot ER distributions
#'
#' Compare the ER score distribution from a collection of new studies to those in the \link{deep_landscape} dataset.
#'
#' @param x \code{\link{deep_mutational_scan}}.
#' @return A \code{\link[ggplot2]{ggplot2}} plot.
#' @examples
#' # Plot a single studies ER
#' dms <- deepscanscape::deep_scans$p53
#' plot_er_distribution(dms)
#'
#' # Plot multiple studies
#' comb_dms <- bind_scans(dms, annotate_missing = TRUE)
#' plot_er_distribution(comb_dms)
#' @export
plot_er_distribution <- function(x) {
  background <- dplyr::select(deepscanscape::deep_landscape, .data$A:.data$Y)
  background <- tidyr::pivot_longer(background, cols = .data$A:.data$Y, names_to = "mut", values_to = "er")
  background$name <- "Background"

  df <- tidyr::pivot_longer(x$data[c("name", amino_acids)], cols = .data$A:.data$Y, names_to = "mut", values_to = "er")

  ggplot2::ggplot() +
    ggplot2::geom_density(data = background, mapping = ggplot2::aes(x = .data$er, fill = "Background"),
                          colour = "grey") +
    ggplot2::stat_density(geom = "line", data = df, position = "identity",
                          mapping = ggplot2::aes(x = .data$er, y = .data$..density.., colour = .data$name)) +
    ggplot2::scale_colour_brewer(name = "", palette = "Set1", direction = "qual") +
    ggplot2::scale_fill_manual(name = "", values = c(Background = "grey")) +
    ggplot2::labs(x = "Normalised ER", y = "Density") +
    theme_deepscanscape()
}

# TODO allow labels for specific points
#' Plot a new study on the deep mutational landscape
#'
#' Project new data onto the deep mutational landscape derived from the combined data, viewed in UMAP space. This view
#' allows you to determine if new data follows the expected distribution and compare it with the average properties for
#' positions in the same regions of the landscape. This can give insights into the likely properties of positions in the
#' new protein.
#'
#' Points represent positions from the new protein(s), positioned based on predicted UMAP coordinates from the
#' combined landscape model. These are plotted on top of the whole combined landscape, which is represented as grey
#' points when no feature is used and coloured hexagonal bins when one is. These bins are coloured based on the mean
#' value of the chosen feature for all points that fall within the hexagon for numerical landscape features.
#' If the chosen feature is a FoldX Gibbs Free Energy prediction values are initially clamped to -5 <= x <= 5 because
#' larger values lead to difficult to read scales and are generally thought to represent program quirks rather than
#' real biology. There are also two special values for feature "count" and "cluster", which colour bins by the
#' number of landscape positions and the most common cluster type, respectively.
#'
#' The input data is annotated using \code{\link{annotate}} if it is not already, based on the \code{annotated} flag
#' for deep_mutational_scan objects and the presence of the required \code{umap1/2} columns for data frames. It is
#' better to save the results of \code{annotate} if doing multiple downstream analyses because can be a relatively slow
#' operation.
#'
#' @param x \code{\link{deep_mutational_scan}}.
#' @param feature String name of a feature to project onto the background landscape. This can be a numeric column from
#' the \code{\link{deep_landscape}} dataset or the special values "count" or "subtype", which map the position density
#' and most common subtype category respectively.
#' @return A \code{\link[ggplot2]{ggplot2}} plot.
#' @examples
#' dms <- deepscanscape::deep_scans$p53
#'
#' # Plot point positions
#' plot_landscape(dms)
#'
#' # Plot against mean conservation
#' dms <- annotate(dms)
#' plot_landscape(dms, feature = "mean_sift")
#'
#' # Plot against spceial features
#' plot_landscape(dms, feature = "count")
#' plot_landscape(dms, feature = "cluster")
#'
#' # Plot multiple studies
#' comb_dms <- bind_scans(dms, annotate_missing = TRUE)
#' plot_landscape(comb_dms, feature = "total_energy")
#'
#' @export
plot_landscape <- function(x, feature = NULL) {
  if (!is.deep_mutational_scan(x)) {
    stop("x is not a deep_mutational_scan")
  }

  if (!x$annotated) {
    warning("deep_mutational_scan is not annotated. Annotating using annotate().", immediate. = TRUE)
    x <- annotate(x)
  }

  if (is.null(feature)) {
    p <- plot_landscape_plain(x)
  } else {
    p <- plot_landscape_feature(x, feature = feature)
  }

  return(p)
}

#' Plot new data onto the plain deep mutational landscape
#'
#' Internal helper called by plot_landscape when no feature is passed
#'
#' @param x \code{\link{deep_mutational_scan}}.
#' @return A \code{\link[ggplot2]{ggplot2}} plot.
#' @keywords internal
plot_landscape_plain <- function(x) {
  # Determine study colours
  n_studies <- nrow(x$meta)
  if (n_studies == 1) {
    col_scale <- ggplot2::scale_colour_manual(name = "", values = "black")
  } else if (n_studies <= 8) {
    col_scale <- ggplot2::scale_colour_brewer(name = "", type = "qual", palette = "Set1")
  } else {
    col_scale <- ggplot2::scale_colour_brewer(name = "", type = "qual", palette = "Paired")
  }

  # Plot
  ggplot2::ggplot(mapping = ggplot2::aes(x = .data$umap1, y = .data$umap2)) +
    ggplot2::geom_point(data = deepscanscape::deep_landscape, mapping = ggplot2::aes(fill = "Background"),
                        shape = 21, colour = "grey") +
    ggplot2::geom_point(data = x$data, mapping = ggplot2::aes(colour = .data$name), shape = 16) +
    ggplot2::scale_fill_manual(name = "", values = c(Background = "grey")) +
    col_scale +
    ggplot2::labs(x = "UMAP1", y = "UMAP2") +
    theme_deepscanscape()
}

#' Plot new data on top of deep mutational landscape features
#'
#' Internal helper called by plot_landscape when called with features.
#'
#' @param x \code{\link{deep_mutational_scan}}.
#' @param feature String name of a feature to project onto the background landscape. This can be a numeric column from
#' the \code{\link{deep_landscape}} dataset or the special values "count" or "cluster", which map the position density
#' and most common subtype category respectively.
#' @return A \code{\link[ggplot2]{ggplot2}} plot.
#' @keywords internal
plot_landscape_feature <- function(x, feature) {
  comb_df <- deepscanscape::deep_landscape

  if (feature %in% c(names(deepscanscape::deep_landscape))) {
    if (typeof(deepscanscape::deep_landscape[[feature]]) == "character" & feature != "cluster") {
      stop("Deep landscape features must be numeric variables (e.g. 'mean_sift' or 'entropy_sidechain') or 'cluster'")
    }

    comb_df <- deepscanscape::deep_landscape[!is.na(deepscanscape::deep_landscape[[feature]]), ]

    if (feature %in% foldx_terms) {
      comb_df[feature] <- clamp(comb_df[[feature]], -5, 5)
    }
  } else if (!feature %in% c("count")) {
    stop("Feature must be a variable in the deep_landscape dataset, 'count' or 'cluster'")
  }

  # Process feature, selecting the fill scale and summary function
  if (feature == "count") {
    func <- length
    fill_scale <- ggplot2::scale_fill_distiller(name = "Positions", type = "seq", palette = "YlGnBu", direction = 1)
    feature <- "mean_sift" # any feature can be used, as all will be the same length
    bins <- 40

  } else if (feature == "cluster") {
    comb_df <- dplyr::left_join(comb_df, deepscanscape::subtypes[c("cluster", "group")], by = "cluster")
    comb_df$group <- as.factor(comb_df$group)
    group_levels <- levels(comb_df$group)
    comb_df$group_num <- as.integer(comb_df$group)

    group_cols <- c(`Aliphatic` = "#ffff33", `Aromatic` = "#377eb8", `Large aliphatic` = "#a65628",
                    `Negative` = "#f781bf", `Not aromatic` = "#e41a1c", `Not proline` = "#4daf4a",
                    `Outlier` = "#000000", `Permissive` = "#999999", `Positive` = "#984ea3",
                    `Small aliphatic` = "#ff7f00", `Unique` = "#00ffff")

    fill_scale <- ggplot2::scale_fill_manual(name = "Subtype Group", values = group_cols)

    func <- function(x) {
      x <- table(x)
      x <- x[names(x) != "7"] # Exclude outliers
      return(group_levels[as.integer(names(x)[which.max(x)])])
    }
    feature <- "group_num"
    bins <- 30
  } else {
    func <- mean
    fill_scale <- get_feature_scale(feature)
    bins <- 40
  }

  # Select colour scale for studies
  n_studies <- nrow(x$meta)
  if (n_studies == 1) {
    col_scale <- ggplot2::scale_colour_manual(name = "", values = "black")
  } else if (n_studies <= 8) {
    col_scale <- ggplot2::scale_colour_brewer(name = "", type = "qual", palette = "Set1")
  } else {
    col_scale <- ggplot2::scale_colour_brewer(name = "", type = "qual", palette = "Paired")
  }

  # Make base plot
  p <- ggplot2::ggplot(mapping = ggplot2::aes(x = .data$umap1, y = .data$umap2)) +
    ggplot2::stat_summary_hex(data = comb_df, fun = func, bins = bins, mapping = ggplot2::aes(z = .data[[feature]])) +
    fill_scale +
    col_scale +
    ggplot2::labs(x = "UMAP1", y = "UMAP2") +
    theme_deepscanscape()

  # Add points depending on shape required
  if (n_studies == 1) {
    p <- p + ggplot2::geom_point(data = x$data, mapping = ggplot2::aes(colour = .data$name), shape = 16)
  } else if (n_studies <= 6) {
    p <- p + ggplot2::geom_point(data = x$data, mapping = ggplot2::aes(shape = .data$name), size = 2.5) +
      ggplot2::geom_point(data = x$data, mapping = ggplot2::aes(shape = .data$name, colour = .data$name), size = 0.9) +
      ggplot2::scale_shape(name = "")
  } else {
    p <- p + ggplot2::geom_point(data = x$data, shape = 16, size = 2.5, colour = "black") +
      ggplot2::geom_point(data = x$data, mapping = ggplot2::aes(colour = .data$name), shape = 16, size = 0.9)
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
    fill_scale <- scale_gradientn(name = "Surface Accessibility", colours = c("#1a2a6c", "#b21f1f", "#fdbb2d"),
                                  values = c(0, 0.4, 1))
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
#' permissive positions a weakly conserved protein. Proportions can also be explicitly compared to the frequencies
#' in the \link{deep_landscape} dataset. In this case differences are tested using an FDR corrected two sided binomial
#' test.
#'
#' @param x \code{\link{deep_mutational_scan}}.
#' @param compare Compare frequencies to those in the base dataset.
#' @return A \code{\link[ggplot2]{ggplot2}} plot. This is a horizontal stacked bar plot when compare = FALSE and a
#' vertical side by side bar plot for compare = TRUE
#' @examples
#' dms <- annotate(deepscanscape::deep_scans$p53)
#' plot_cluster_frequencies(dms)
#'
#' # Plot multiple studies
#' comb_dms <- bind_scans(dms, annotate_missing = TRUE)
#' plot_cluster_frequencies(comb_dms)
#'
#' # Compare to the deep_landscape dataset
#' plot_cluster_frequencies(comb_dms, compare = TRUE)
#'
#' @export
plot_cluster_frequencies <- function(x, compare = FALSE) {
  if (!is.deep_mutational_scan(x)) {
    stop("x is not a deep_mutational_scan")
  }

  if (!x$annotated) {
    warning("deep_mutational_scan is not annotated. Annotating using annotate().", immediate. = TRUE)
    x <- annotate(x)
  }

  new_freq <- dplyr::summarise(dplyr::group_by(x$data, .data$wt, .data$cluster), n = dplyr::n(), .groups = "drop_last")
  new_freq <- dplyr::ungroup(dplyr::mutate(new_freq, prop = .data$n / sum(.data$n)))
  new_freq$wt <- factor(new_freq$wt, levels = sort(unique(new_freq$wt), decreasing = TRUE))
  new_freq$cluster <- stringr::str_sub(new_freq$cluster, start = 2)

  if (compare) {
    background <- dplyr::summarise(dplyr::group_by(deepscanscape::deep_landscape, .data$wt, .data$cluster),
                                   n = dplyr::n(), .groups = "drop_last")
    background <- dplyr::ungroup(dplyr::mutate(background, prop = .data$n / sum(.data$n)))
    background$wt <- factor(background$wt, levels = sort(unique(new_freq$wt), decreasing = TRUE))
    background$cluster <- stringr::str_sub(background$cluster, start = 2)

    overview <- dplyr::bind_rows(Background = background, `New Data` = new_freq, .id = "type")
    overview$cluster <- factor(overview$cluster, levels = c(seq_len(8), "P", "O", "A"))
    overview <- tidyr::complete(overview, .data$wt, .data$cluster, .data$type, fill = list(n = 0, prop = 0))

    # Test new freqs
    binom <- new_freq[new_freq$cluster != "A", c("wt", "cluster", "n", "prop")]
    binom <- dplyr::group_by(binom, .data$wt)
    binom <- dplyr::ungroup(dplyr::mutate(binom, tot = sum(.data$n)))
    binom <- dplyr::left_join(binom, dplyr::select(background, .data$wt, .data$cluster, p = .data$prop),
                              by = c("wt", "cluster"))
    binom$pvalue <- mapply(function(x, n, p) stats::binom.test(x, n, p)$p.value,
                           x = binom$n, n = binom$tot, p = binom$p)
    binom$padj <- stats::p.adjust(binom$pvalue, method = "fdr")
    binom$symb <- as.character(stats::symnum(binom$padj, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                                             symbols = c("***", "**", "*", ".", " ")))

    p <- ggplot2::ggplot(overview, ggplot2::aes(x = .data$cluster, y = .data$prop)) +
      ggplot2::facet_wrap(~wt, nrow = 4, strip.position = "bottom", scales = "free_x") +
      ggplot2::scale_y_continuous(labels = function(x) stringr::str_c(100 * x, "%"), breaks = seq(0, 1, 0.25),
                                  limits = c(0, 1.15)) +
      ggplot2::geom_col(ggplot2::aes(fill = .data$type), position = "dodge") +
      ggplot2::geom_text(data = binom, mapping = ggplot2::aes(label = .data$symb, y = .data$prop + 0.15)) +
      ggplot2::scale_fill_manual(name = "", values = c(Background = "black", `New Data` = "red")) +
      ggplot2::labs(x = "", y = "Percentage of Positions") +
      ggplot2::coord_cartesian(clip = "off") +
      theme_deepscanscape() +
      ggplot2::theme(strip.placement = "outside")

  } else {
    max_cluster <- max(as.integer(stringr::str_subset(new_freq$cluster, "[0-9]+")))
    new_freq$cluster <- factor(new_freq$cluster, levels = c("O", "A", "P", rev(seq_len(max_cluster))))

    counts <- dplyr::summarise(dplyr::group_by(new_freq, .data$wt), n = sum(.data$n), .groups = "drop")
    sec_axis <- ggplot2::dup_axis(name = "", labels = counts$n)

    if (max_cluster <= 8) {
      fill_scale <- ggplot2::scale_fill_manual(values = c("1" = "#e41a1c", "2" = "#377eb8", "3" = "#4daf4a",
                                                          "4" = "#984ea3", "5" = "#ff7f00", "6" = "#ffff33",
                                                          "7" = "#42b7ce", "8" = "#f781bf",
                                                          "P" = "#adadad", "A" = "#a65628", "O" = "#666666"),
                                               name = "Subtype")
    } else {
      fill_scale <- ggplot2::scale_fill_brewer(name = "Subtype")
    }

    p <- ggplot2::ggplot(new_freq, ggplot2::aes(x = as.integer(.data$wt), y = .data$prop, fill = .data$cluster)) +
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
  return(p)
}
