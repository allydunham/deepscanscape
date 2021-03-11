# Recluster the whole dataset with new data and report on the output results

# TODO add more parameters to clustering?
#' Recluster DMS data with a combined dataset
#'
#' @param x data frame or list of \link{deep_mutational_scan} objects to recluster
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
