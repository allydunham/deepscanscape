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
#' @param x \link{deep_mutational_scan} to analyse
#' @export
plot_dms_heatmap <- function(x) {
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

#' Plot a new study on the deep mutational landscape
#'
#' @param x \link{deep_mutational_scan} to analyse. Unannotated datasets will be annotated using \link{annotate_dms},
#'   but it is generally better to do this beforehand as it an expensive operation.
#' @param name Name of this dataset. When NULL a name is generated from the object data
#' @export
plot_dms_landscape <- function(x, name = NULL) {
  if (!x$annotated) {
    warning("deep_mutational_scan is not annotated. Annotating using annotate_dms().")
    x <- annotate_dms(x)
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

  ggplot2::ggplot(mapping = ggplot2::aes(x = .data$umap1, y = .data$umap2)) +
    ggplot2::geom_point(data = deepscanscape::deep_mutational_scans, mapping = ggplot2::aes(colour = "background")) +
    ggplot2::geom_point(data = x$data, mapping = ggplot2::aes(colour = "this")) +
    ggplot2::scale_colour_manual(values = c(background = "grey", this = "red"), name = "",
                                 labels = c(background = "Background", this = name)) +
    ggplot2::labs(x = "UMAP1", y = "UMAP2") +
    theme_deepscanscape()
}
