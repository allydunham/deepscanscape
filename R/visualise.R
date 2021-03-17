# Visualise new DMS datasets
# TODO - test

#' Default deepscanscape GGPlot2 theme
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

#' Plot a deep mutational scan dataset as a heatmap
#'
#' @param x \code{\link{deep_mutational_scan}} to analyse
#' @export
plot_dms_heatmap <- function(x) {
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

# TODO Update to support combined data
#' Plot a new study on the deep mutational landscape
#'
#' @param x \link{deep_mutational_scan} to analyse. Unannotated datasets will be annotated using \link{annotate},
#'   but it is generally better to do this beforehand as it an expensive operation.
#' @param name Name of this dataset. When NULL a name is generated from the object data
#' @param feature Feature from the \code{\link{deep_mutational_scan}} dataset to map positions against
#' @export
plot_landscape <- function(x, name = NULL, feature = NULL) {
  if (!x$annotated) {
    warning("deep_mutational_scan is not annotated. Annotating using annotate().")
    x <- annotate(x)
  }

  # Guess good
  if (is.null(name)) {
    gene_na <- is.na(x$gene) | is.null(x$gene)
    study_na <- is.na(x$study) | is.null(x$study)
    if (gene_na & study_na) {
      name <- "This study"
    } else if (study_na) {
      name <- x$gene
    } else {
      name <- x$study
    }
  }

  if (is.null(feature)) {
    p <- ggplot2::ggplot(mapping = ggplot2::aes(x = .data$umap1, y = .data$umap2)) +
      ggplot2::geom_point(data = deepscanscape::deep_mutational_scans, mapping = ggplot2::aes(colour = "background")) +
      ggplot2::geom_point(data = x$data, mapping = ggplot2::aes(colour = "this")) +
      ggplot2::scale_colour_manual(values = c(background = "grey", this = "red"), name = "",
                                   labels = c(background = "Background", this = name)) +
      ggplot2::labs(x = "UMAP1", y = "UMAP2") +
      theme_deepscanscape()
  } else {
    if (!feature %in% names(deepscanscape::deep_mutational_scans)) {
      stop("Feature must be a variable in the deep_mutational_scans dataset")
    }

    if (typeof(deepscanscape::deep_mutational_scans[[feature]]) == "character") {
      stop("Features must be numeric variables (e.g. 'mean_sift' or 'entropy_sidechain')")
    }

    pretty <- stringr::str_to_title(stringr::str_replace_all(feature, "_", " "))
    df <- deepscanscape::deep_mutational_scans[!is.na(deepscanscape::deep_mutational_scans[[feature]]), ]
    if (feature == "mean_sift") {
      fill_scale <- ggplot2::scale_fill_distiller(name = expression(log[10] * "SIFT"), type = "seq", palette = "RdPu")
    } else if (feature == "mean_score") {
      fill_scale <- ggplot2::scale_fill_gradient2(name = "Mean ER", low = "#d73027", mid = "#ffffbf", high = "#4575b4")
    } else if (feature %in% amino_acids) {
      fill_scale <- ggplot2::scale_fill_gradient2(name = feature, low = "#4575b4", mid = "#ffffbf", high = "#d73027")
    } else if (feature == "all_atom_rel") {
      fill_scale <- ggplot2::scale_fill_gradientn(colours = c("#1a2a6c", "#b21f1f", "#fdbb2d"), values = c(0, 0.2, 1))
    } else if (feature %in% c("total_energy", "backbone_hbond", "sidechain_hbond", "van_der_waals", "electrostatics",
                              "solvation_polar", "solvation_hydrophobic", "van_der_waals_clashes", "entropy_sidechain",
                              "entropy_mainchain", "cis_bond", "torsional_clash", "backbone_clash", "helix_dipole",
                              "disulfide", "electrostatic_kon", "partial_covalent_bonds", "energy_ionisation")) {
      # TODO document clamping
      df[feature] <- clamp(df[[feature]], -5, 5)
      fill_scale <- ggplot2::scale_fill_gradient2(name = pretty, low = "#4575b4",
                                                  mid = "#ffffbf", high = "#d73027")
    } else if (feature == "hydrophobicity") {
      fill_scale <- ggplot2::scale_fill_gradientn(colours = c("#4575b4", "#e0f3f8", "#fee090", "#fc8d59", "#d73027"),
                                                    values = scales::rescale(c(-0.4, 0, 0.4, 0.8, 1.2)),
                                                    limits = c(-0.4, 1.201))
    } else if (feature %in% c("phi",  "psi")) {
      fill_scale <- ggplot2::scale_fill_gradient2(name = pretty, low = "#c51b7d", mid = "#ffffbf", high = "#4d9221")
    } else {
      fill_scale <- ggplot2::scale_fill_continuous(name = pretty)
    }

    p <- ggplot2::ggplot(mapping = ggplot2::aes(x = .data$umap1, y = .data$umap2)) +
      ggplot2::stat_summary_hex(data = df, fun = mean, bins = 40, mapping = ggplot2::aes(z = .data[[feature]])) +
      ggplot2::geom_point(data = x$data, mapping = ggplot2::aes(colour = "this")) +
      ggplot2::scale_colour_manual(name = "", values = c(this = "black"), labels = c(this = name)) +
      fill_scale +
      ggplot2::labs(x = "UMAP1", y = "UMAP2") +
      theme_deepscanscape()
  }
  return(p)
}

#' Plot amino acid subtype frequencies
#'
#' @param x \code{\link{deep_mutational_scan}} or \link[=rbind.deep_mutational_scan]{combined DMS data frame}.
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
  overview$cluster_num <- factor(overview$cluster_num, levels = c("O", "P", seq_len(max_cluster)))

  counts <- dplyr::summarise(dplyr::group_by(overview, .data$wt), n = sum(.data$n), .groups = "drop")
  sec_axis <- ggplot2::dup_axis(name = "", labels = counts$n)

  if (max_cluster <= 7) {
    pal <- "Set1"
  } else {
    if (max_cluster > 10) {
      warning(">10 clusters for at least one amino acid, colours may not be ideal")
    }
    pal <- "Set3"
  }

  ggplot2::ggplot(overview, ggplot2::aes(x = as.integer(.data$wt), y = .data$prop, fill = .data$cluster_num)) +
    ggplot2::scale_x_continuous(breaks = seq_len(20), labels = counts$wt, sec.axis = sec_axis) +
    ggplot2::scale_y_continuous(labels = function(x) stringr::str_c(100 * x, "%")) +
    ggplot2::coord_flip(expand = FALSE) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::labs(x = "WT Amino Acid", y = "Percentage of Positions") +
    ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE)) +
    ggplot2::scale_fill_brewer(name = "Subtype", type = "qual", palette = pal) +
    theme_deepscanscape() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank())
}
