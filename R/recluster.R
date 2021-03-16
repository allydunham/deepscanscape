# Recluster the whole dataset with new data and report on the output results

# TODO add more parameters to clustering?
#' Recluster DMS data with a combined dataset
#'
#' @param x data frame or list of \code{\link{deep_mutational_scan}} objects to recluster
#' @param keep_clustering logical. Keep the outputs of \link[stats]{hclust} and \link[dynamicTreeCut]{cutreeHybrid} for
#' downstream analysis
#' @param deep_split Named vector of deepSplit parameters to pass to \link[dynamicTreeCut]{cutreeHybrid}. Must be a
#' numeric vector with a named entry for each amino acid or a single integer to apply to all amino acids
#' @param add_combined Combine the supplied data with the data contained in \link{deep_mutational_scans}
#'
#' @export
recluster_dms <- function(x, keep_clustering = FALSE, deep_split=NULL, add_combined=TRUE) {
  # Process deepSplit
  if (is.null(deep_split)) {
    deep_split <- c("A" =  0, "C" =  0, "D" =  1, "E" =  0, "F" =  1, "G" =  1, "H" =  0, "I" =  0, "K" =  1,
                    "L" =  0, "M" =  0, "N" =  0, "P" =  1, "Q" =  1, "R" =  1, "S" =  0, "T" =  1, "V" =  0,
                    "W" =  0, "Y" =  1)
  } else {
    if (!all(deep_split %in% 0:4)) {
      stop("Invalid deepSplit value. Must be an integer or a named vector with values for each amino acid.",
           "Each value must be a integer from 0 to 4")
    }
    if (length(deep_split) == 1) {
      deep_split <- structure(rep(deep_split, 20), names = amino_acids)
    } else {
      if (!(length(deep_split) == 20 & all(sort(names(deep_split)) == amino_acids))) {
        stop("Invalid deepSplit value. Must be an integer or a named vector with values for each amino acid.",
             "Each value must be a integer from 0 to 4")
      }
    }
  }

  # Generate data
  if (is.deep_mutational_scan(x)) {
    df <- tibble::as_tibble(x, full = TRUE)
  } else if (is.data.frame(x)) {
    if (!all(c("study", "gene", "position", "wt", amino_acids) %in% names(x))) {
      stop("x does not contain the required columns (see ?recluster_dms)")
    }
    df <- tibble::as_tibble(x)
  } else if (inherits(x, "list")) {
    if (!all(sapply(x, is.deep_mutational_scan))) {
      stop("Datasets to recluster must all be deep mutational scans.\nGenerate these using deep_mutational_scan()")
    }
    df <- do.call(rbind, x)
  } else {
    stop("Unrecognised input\n Pass a correctly formatted data frame or a list of deep_mutational_scan() objects")
  }

  df <- df[c("study", "gene", "position", "wt", amino_acids)]
  if (add_combined) {
    comb <- dplyr::select(deepscanscape::deep_mutational_scans, .data$study, .data$gene,
                          .data$position, .data$wt, dplyr::one_of(amino_acids))
    df <- dplyr::bind_rows(df, comb)
  }

  # Perform PCA
  er_mat <- as.matrix(df[amino_acids])
  pca <- stats::prcomp(er_mat)
  df <- dplyr::bind_cols(df, tibble::as_tibble(pca$x))

  # Identify Permissive Positions
  permissive <- apply(abs(er_mat) < 0.4, 1, all)
  df_permissive <- df[permissive, ]
  df_permissive$cluster <- stringr::str_c(df_permissive$wt, "P")

  # Clustering
  df <- df[!permissive, ]
  aa_dfs <- split(df, f = df$wt, drop = TRUE)

  cluster_wrapper <- function(aa) {
    return(cluster_df(aa_dfs[[aa]], deep_split = deep_split[aa]))
  }
  clusts <- sapply(names(aa_dfs), cluster_wrapper, simplify = FALSE)

  # Assemble output
  df <- dplyr::bind_rows(lapply(clusts, "[[", "tbl"))
  df$cluster <- stringr::str_c(df$wt, df$cluster)
  df <- dplyr::bind_rows(df, df_permissive) # Done in two stages so R recognises one is a list of D
  df <- dplyr::select(df, .data$cluster, dplyr::everything())
  df <- dplyr::arrange(df, .data$study, .data$gene, .data$position)

  # Rename outlier clusters
  outlier_inds <- stringr::str_detect(df$cluster, "[A-Z]0")
  df$cluster[outlier_inds] <- stringr::str_replace(df$cluster[outlier_inds], "0", "O")

  if (keep_clustering) {
    return(list(data = df, hclust = lapply(clusts, "[[", "hclust"), treeCut = lapply(clusts, "[[", "treeCut")))
  }
  return(df)
}

# TODO Allow varying columns?
#' Cluster rows of a data frame using hierarchical clustering and dynamic tree cutting
#'
#' @param tbl Data frame containing columns
#' @param deep_split deepSplit parameter passed to \link[dynamicTreeCut]{cutreeHybrid}
cluster_df <- function(tbl, deep_split = 0) {
  m <- as.matrix(tbl[stringr::str_c("PC", 2:20)])
  d <- stats::as.dist(cosine_distance_matrix(m))
  hc <- stats::hclust(d = d, method = "average")
  invisible(utils::capture.output(
    clus <- dynamicTreeCut::cutreeHybrid(dendro = hc, distM = as.matrix(d), deepSplit = deep_split)
  ))
  tbl$cluster <- clus$labels
  return(list(tbl = tbl, hclust = hc, treeCut = clus))
}

# TODO export these functions? maybe useful if people want to tinker with clustering
#' Cosine distance matrix
#'
#' Calculate the cosine distance between all rows of a matrix.
#'
#' @param m Numeric matrix
cosine_distance_matrix <- function(m) {
  # TODO - document this properly. Rounding is required to account for floating point precision
  cosine <- tcrossprod(m) / sqrt(tcrossprod(rowSums(m^2)))
  cosine <- acos(round(cosine, digits = 8)) / pi
  return(cosine)
}

#' Summarise a deep mutational scan recluster
#'
#' @param df reclustered data, as output by \link{recluster_dms}
#' @param aa Character vector of amino acids to summarise clusters from.
#' @param square_tiles Force heatmap tiles to be square. It can be useful to disable this when summarising a large
#'   number of clusters.
plot_dms_recluster <- function(df, aa = NULL, square_tiles = TRUE) {
  req_cols <- c("cluster", "study",  "gene", "position", "wt", amino_acids)
  if (!all(req_cols %in% names(df))) {
    stop("Incorrect columns in df, must include: ", stringr::str_c(req_cols, collapse = ", "))
  }

  if (!is.null(aa)) {
    if (!all(aa %in% amino_acids)) {
      stop("Incorrect AA argument")
    }
    df <- df[df$wt %in% aa, ]
  }

  # Calculate cluster profiles and counts
  profile <- dplyr::summarise(dplyr::group_by(df, .data$cluster), dplyr::across(dplyr::one_of(amino_acids), mean))
  profile <- tidyr::pivot_longer(profile, -.data$cluster, names_to = "mut", values_to = "er")
  profile$mut <- as.factor(profile$mut)

  counts <- dplyr::summarise(dplyr::group_by(df, .data$cluster), n = dplyr::n())
  counts$percent <- counts$n / sum(counts$n)

  # Set cluster order
  cluster_levels <- unique(profile$cluster)
  cluster_levels <- c(sort(stringr::str_subset(cluster_levels, "[A-Z]O"), decreasing = TRUE),
                      sort(stringr::str_subset(cluster_levels, "[A-Z]P"), decreasing = TRUE),
                      sort(stringr::str_subset(cluster_levels, "[A-Z][0-9]+"), decreasing = TRUE))

  profile$cluster <- factor(profile$cluster, levels = cluster_levels)
  counts$cluster <- factor(counts$cluster, levels = levels(profile$cluster))

  # Format axes
  y_labels <- levels(profile$cluster)
  y_sec_labels <- counts$n[order(counts$cluster)]

  if (max(counts$percent) < 0.05) {
    x_labels <- c(amino_acids, "", "0%", "1.25%", "2.5%", "3.75%", "5%")
    x_mult <- 80
  } else if (max(counts$percent) < 0.1) {
    x_labels <- c(amino_acids, "", "0%", "2.5%", "5%", "7.5%", "10%")
    x_mult <- 40
  } else if (max(counts$percent) < 0.2) {
    x_labels <- c(amino_acids, "", "0%", "5%", "10%", "15%", "20%")
    x_mult <- 20
  } else if (max(counts$percent) < 0.4) {
    x_labels <- c(amino_acids, "", "0%", "10%", "20%", "30%", "40%")
    x_mult <- 10
  } else {
    x_labels <- c(amino_acids, "", "0%", "25%", "50%", "75%", "100%")
    x_mult <- 4
  }

  # Calculate aesthetic layers
  tile_aes <- ggplot2::aes(x = as.integer(.data$mut), y = as.integer(.data$cluster), fill = .data$er)
  bar_aes <- ggplot2::aes(xmin = 22, xmax = 22 + x_mult * .data$percent, ymin = as.integer(.data$cluster) - 0.25,
                          ymax = as.integer(.data$cluster) + 0.25)

  # Make plot
  p <- ggplot2::ggplot() +
    ggplot2::scale_x_continuous(breaks = 1:26, limits = c(0, 26), labels = x_labels) +
    ggplot2::scale_y_continuous(breaks = seq_along(y_labels), limits = c(0.5, length(y_labels) + 0.5),
                                labels = y_labels, sec.axis = ggplot2::dup_axis(name = "", labels = y_sec_labels)) +
    ggplot2::geom_vline(xintercept = 22) +
    ggplot2::geom_vline(xintercept = 23:26, linetype = "dotted") +
    ggplot2::geom_tile(data = profile, mapping = tile_aes) +
    ggplot2::geom_rect(data = counts, mapping = bar_aes, fill = "#542788") +
    ggplot2::scale_fill_gradient2(name = "Mean ER", low = "#a50026", mid = "#ffffbf", high = "#313695") +
    theme_deepscanscape() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   axis.title = ggplot2::element_blank())

  if (square_tiles) {
    p <- p + ggplot2::coord_equal(expand = FALSE)
  } else {
    p <- p + ggplot2::coord_cartesian(expand = FALSE)
  }
  return(p)
}
