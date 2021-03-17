# S3 class to store a dms dataset
# TODO - Test these funcs

#' Internal deep_mutational_scan constructor
#'
#' Create a new raw deep_mutational_scan object. This is expected to have no annotations and no imputation but have
#' had scores transformed and normalised onto the standard scale.
#'
#' @param df A data frame in the required format (position, wt, A-Y, additional data...).
#' @param study The source study.
#' @param gene The studied gene.
#' @return A deep_mutational_scan S3 object (see \code{\link{deep_mutational_scan}})
new_deep_mutational_scan <- function(df, study, gene) {
  out <- list(data = tibble::as_tibble(df), study = study, gene = gene,
              annotated = FALSE, impute_mask = NA, cluster = NA)
  class(out) <- c("deep_mutational_scan")
  return(out)
}

# TODO Use this check more widely?
# TODO Export this?
#' Validate deep_mutational_scan objects
#'
#' Check properties of deep_mutational_scan object adhere to various baseline expectations, which would lead to
#' downstream errors if not caught.
#'
#' The following assumptions are currently checked:
#' \itemize{
#'   \item Correct class set
#'   \item Correct fields present
#'   \item \code{data} is a tibble
#'   \item \code{annotated} is logical
#'   \item \code{impute_matrix} is a tibble or NA
#'   \item \code{cluster} is a tibble or NA
#'   \item \code{data} contains the correct columns
#'   \item \code{data} contains no duplicate rows
#'   \item \code{impute_matrix} contains the correct type of values
#'   \item Annotated objects contain the correct \code{data} and \code{cluster} columns
#' }
#'
#' @param x \code{\link{deep_mutational_scan}} object.
#' @return The input is returned unaltered, since the main purpose of the function is to raise an error if the
#'   object is erroneous.
validate_deep_mutational_scan <- function(x) {
  # Check fundamental structure
  if (!is.deep_mutational_scan(x)) {
    stop("Not a deep_mutational_scan (check class(x))")
  }

  if (!all(c("data", "study", "gene", "annotated", "impute_mask", "cluster") %in% names(x))) {
    stop("Missing field(s). The object must contain 'data', 'study', 'gene', 'annotated' & 'impute_mask'")
  }

  # Check field types
  if (!tibble::is_tibble(x$data)) {
    stop("Invalid 'data'. Must be a tibble")
  }

  if (!is.logical(x$annotated)) {
    stop("Invalid 'annotated' value (", x$annotated, "). Must be TRUE/FALSE")
  }

  if (!((is.matrix(x$impute_mask) & is.numeric(x$impute_mask)) | identical(x$impute_mask, NA))) {
    stop("Invalid 'impute_mask'. Must be a logical matrix or NA")
  }

  if (!(identical(x$cluster, NA) | tibble::is_tibble(x$cluster))) {
    stop("Invalid 'cluster'. Must be a tibble or NA")
  }

  # Check data object
  missing_headers <- c("position", "wt", amino_acids)[which(!c("position", "wt", amino_acids) %in% names(x$data))]
  if (length(missing_headers) > 0) {
    stop("Missing data columns: ", paste(headers, sep = ", "))
  }

  if (nrow(x$data) != nrow(dplyr::distinct(x$data[c("position", "wt")]))) {
    stop("Data contains duplicated positions")
  }

  # TODO - Check scores fit right profile? (Maybe too fragile)

  # Check impute_matrix
  if (!identical(x$impute_mask, NA)) {
    if (!all(unique(c(x$impute_mask)) %in% 0:2)) {
      stop("Impute matrix contains incorrect values: must all be 0, 1 or 2")
    }
  }

  # Check annotation
  if (x$annotated) {
    headers <- c("cluster", stringr::str_c("PC", seq_len(20)), "umap1", "umap2")
    missing_headers <- headers[which(!headers %in% names(x$data))]
    if (length(missing_headers) > 0) {
      stop("Missing annotation columns: ", paste(missing_headers, sep = ", "))
    }

    if (identical(x$cluster, NA)) {
      stop("Cluster annotation missing")
    }

    headers <- c("cluster", "base_cluster", "permissive", "small_margin", "high_distance",
                 "dist1", "dist2", "dist3", "dist4", "dist5", "dist6", "dist7", "dist8", "notes")
    missing_headers <- headers[which(!headers %in% names(x$cluster))]
    if (length(missing_headers) > 0) {
      stop("Missing cluster annotation columns: ", paste(missing_headers, sep = ", "))
    }
  }

  return(x)
}

#' Deep mutational scan data
#'
#' Store deep mutational scanning data in a standardised format, alongside
#' metadata describing the state of the data.
#'
#' @param df A data frame to parse.
#' @param scheme Original data scheme (see \code{\link{parse_dms_data}}).
#' @param trans Function to transform scores onto the standard scale. Accepts a string corresponding to known transforms
#' or a custom function (\code{\link{transform_dms}}).
#' @param na_value Value to set missense NA scores to (see \code{\link{impute_dms}}).
#' @param annotate Annotate the dataset with mutational landscape data (PCA, UMAP and amino acid subtypes).
#' @param study Source of the deep mutational scan.
#' @param gene Gene scanned.
#' @return A deep_mutational_scan S3 object, which broadly behaves as a list containing the following fields:
#'   \itemize{
#'     \item data: \code{\link[tibble]{tibble}} containing positions ER scores and other positional data
#'     \item study: source study
#'     \item gene: gene scanned
#'     \item annotated: logical indicating if the data has been annotated with PCs, UMAP coordinates and clusters
#'     \item impute_mask: numeric matrix indicating imputed ER scores, corresponding to the
#'     rows/ER columns of \code{data} (see \code{\link{impute_dms}} for details)
#'     \item cluster: \code{\link[tibble]{tibble}} containing further details on how clusters were assigned during
#'     annotation, corresponding to the rows of \code{data}
#'   }
#'   The class is used like a list apart from [ accesses the main data frame (see \link{dms_extract} for details).
#' @examples
#' path <- system.file("extdata", "urn_mavedb_00000036-a-1_scores.csv",
#'                     package = "deepscanscape")
#' csv <- read.csv(path, skip = 4)
#' dms <- deep_mutational_scan(csv, scheme = "mave", trans = "vamp",
#'                             gene = "LDLRAP1", study = "urn:mavedb:00000036")
#' @export
deep_mutational_scan <- function(df, scheme=NULL, trans=NULL, na_value="impute",
                                 annotate=TRUE, study=NA, gene=NA) {
  # Parse scheme
  if (!is.null(scheme)) {
    df <- parse_dms_data(df, scheme)
  } else {
    df <- tibble::as_tibble(df)
    if (!all(c("position", "wt", "mut", "score") %in% names(df))) {
      stop("Incorrect input df, must include position, wt, mut and score columns")
    }
    df <- dplyr::select(df, .data$position, .data$wt, .data$mut, .data$score, dplyr::everything())
  }

  # Transform
  if (!is.null(trans)) {
    df$score <- transform_dms(df$score, trans)
  }

  # Normalise
  df$score <- normalise_dms(df$score)

  # Format tibble
  df <- df[df$mut %in% amino_acids, ]
  df <- dplyr::arrange(df, .data$position, .data$mut)
  df <- tidyr::pivot_wider(df, names_from = "mut", values_from = "score")
  df <- dplyr::select(df, .data$position, .data$wt, dplyr::all_of(amino_acids), dplyr::everything())

  # Construct object
  out <- validate_deep_mutational_scan(new_deep_mutational_scan(df = df, study = study, gene = gene))

  # Impute
  # TODO note what happens when not doing this
  if (!(is.na(na_value) | is.null(na_value))) {
    out <- impute_dms(out)
  }

  # Annotate
  if (annotate) {
    out <- annotate_dms(out)
  }

  return(validate_deep_mutational_scan(out))
}

#' Determine if an object is a deep_mutational_scan
#'
#' @param x Object to check
#'
#' @export
is.deep_mutational_scan <- function(x) { # nolint
  return(inherits(x, "deep_mutational_scan"))
}

# TODO Document this fully - [ and [[ go to data tbl, $ goes to list
#' Extracting and replacing deep mutational scanning data
#'
#' Extracting and replacing data from \code{\link{deep_mutational_scan}} objects uses a mixture of list and data frame
#' syntax. $ and [[ extract values from the main data fields (e.g. the gene or where the data is annotated) and [
#' provides a shortcut to access the main data table, which can otherwise be accessed via \code{x$data}.
#'
#' @param x \code{\link{deep_mutational_scan}} object.
#' @param i,j,... Indices to access.
#' @param drop Coerce result to lowest possible dimension.
#' @param value Value to set.
#' @examples
#' # Setup object
#' path <- system.file("extdata", "urn_mavedb_00000036-a-1_scores.csv",
#'                     package = "deepscanscape")
#' csv <- read.csv(path, skip = 4)
#' dms <- deep_mutational_scan(csv, scheme = "mave", trans = "vamp",
#'                             gene = "LDLRAP1", study = "urn:mavedb:00000036")
#'
#' # Extract meta data:
#' dms$study
#' dms$impute_mask
#'
#' dms[["gene"]]
#'
#' # Replace meta data
#' dms$gene <- "new_gene"
#' dms$impute_mask <- NA
#'
#' # Quickly access data columns:
#' dms["A"]
#' dms[c(1,2,3), "C"]
#' dms[c("position", "wt", "A", "D")]
#'
#' # Quickly modify data columns
#' dms["position"] <- dms["position"] + 1
#'
#' @name dms_extract
NULL
#> NULL

#' @describeIn dms_extract Extract
#' @export
`[.deep_mutational_scan` <- function(x, i, j, drop = FALSE, ...) {  # nolint
  if (missing(j)) {
    return(x$data[i])
  }
  return(x$data[i, j, ..., drop = drop])
}

#' @describeIn dms_extract Assign
#' @export
`[<-.deep_mutational_scan` <- function(x, i, j, ..., value) {  # nolint
  if (missing(j)) {
    x$data[i] <- value
  } else {
    x$data[i, j, ...] <- value
  }
  return(x)
}

#' General S3 Methods for Deep Mutational Scans
#'
#' \code{\link{deep_mutational_scan}} S3 objects support most generics that make sense for the data structure. Mostly
#' these operate on the underlying list object, although specific \code{format}, \code{print}, \code{str} and \code{dim}
#' methods are provides to make working with the object more convenient. The \code{dim} method provides a direct
#' short-cut to the main data frame.
#'
#' @param x,object \code{\link{deep_mutational_scan}} object.
#' @param indent.str Indent string for \link[utils]{str}.
#' @param nest.lev Nest level for \link[utils]{str}.
#' @param ... Additional arguments.
#'
#' @name dms_s3
NULL
#> NULL

# TODO - Pretty print when extra columns added
# TODO - too many columns?
#' @describeIn dms_s3 S3 format method
#' @export
format.deep_mutational_scan <- function(x, ...) {
  out <- c(paste("# A deep_mutational_scan"),
           paste("# Study:", x$study),
           paste("# Gene:", x$gene),
           paste("#", nrow(x$data), "positions"))

  if (x$annotated) {
    out <- c(out, "# Annotated with clusters, pricipal components and UMAP coordinates")
  }

  out <- c(out, "# Positional data:")

  if (requireNamespace("crayon", quietly = TRUE)) {
    out <- crayon::make_style("darkgrey")(out)
  }

  tbl <- x$data
  tbl <- dplyr::select(tbl, dplyr::any_of("cluster"), .data$position, .data$wt,
                       dplyr::all_of(amino_acids), dplyr::everything())
  tbl_str <- format(tbl)
  out <- c(out, tbl_str[-1])

  return(out)
}

#' @describeIn dms_s3 S3 print method
#' @export
print.deep_mutational_scan <- function(x, ...) { # nolint
  cat(format(x, ...), sep = "\n")
}

#' @describeIn dms_s3 S3 str method
#' @export
str.deep_mutational_scan <- function(object, ..., indent.str = " ", nest.lev = 0) { # nolint
  if (nest.lev != 0L) {
    cat(indent.str)
  }

  cat("deep_mutational_scan [", nrow(object$data), " positions]", "\n", sep = "")

  utils::str(list(study = object$study, gene = object$gene, annotated = object$annotated, data = as.list(object$data)),
             no.list = TRUE, ..., nest.lev = nest.lev + 2L, indent.str = indent.str)
}

#' @describeIn dms_s3 S3 dim method
#' @export
dim.deep_mutational_scan <- function(x) {
  return(dim(x$data))
}

#' Convert deep_mutational_scan objects
#'
#' Conveniently convert \code{\link{deep_mutational_scan}} objects to a \code{\link[tibble]{tibble}} or
#' \code{\link[base]{data.frame}} or \code{\link[base]{list}}, with the option to add gene and study columns.
#'
#' @param x \code{\link{deep_mutational_scan}} object.
#' @param full Include columns containing metadata (study, gene).
#' @param ... Ignored.
#' @return A \code{\link[tibble]{tibble}}, \code{\link[base]{data.frame}} or \code{\link[base]{list}}
#' @importFrom tibble as_tibble
#' @name as_tibble
#' @export
as_tibble.deep_mutational_scan <- function(x, ..., full=FALSE) { # nolint
  if (full) {
    out <- tibble::as_tibble(x)
    out$study <- x$study
    out$gene <- x$gene
    return(dplyr::select(out, .data$study, .data$gene, dplyr::everything()))
  }
  return(x$data)
}

#' @describeIn as_tibble S3 as.data.frame method
#' @export
as.data.frame.deep_mutational_scan <- function(x, ..., full=FALSE) {
  return(as.data.frame(as_tibble(x, full = full, ...)))
}

#' @describeIn as_tibble S3 as.list method
#' @export
as.list.deep_mutational_scan <- function(x, ...) {
  class(x) <- "list"
  return(x)
}

#' Summary plot for Deep Mutational Scans
#'
#' Produce a summary plot for a \code{\link{deep_mutational_scan}} object, showing a heatmap with the ER score for
#' each substitution at each position.
#'
#' @param object,x \code{\link{deep_mutational_scan}} to produce a summary plot for.
#' @param ... Additional parameters.
#' @return A \code{\link[ggplot2]{ggplot2}} object for autoplot or nothing for plot, which prints the object directly.
#'
#' @export
#' @name dms_plot
#' @importFrom ggplot2 autoplot
autoplot.deep_mutational_scan <- function(object, ...){ # nolint
  plot_dms_heatmap(object)
}

#' @describeIn dms_plot S3 plot method
#' @export
#' @importFrom graphics plot
plot.deep_mutational_scan <- function(x, ...) {
  print(autoplot(x, ...))
}

# TODO Move the combined stuff to own file
# TODO support these in all applicable functions
# TODO support passing a list of scans
#' Combine deep mutational scan data
#'
#' Combine multiple \code{\link{deep_mutational_scan}} objects into a single \code{\link[tibble]{tibble}}, including
#' columns giving the study and gene of each position.
#'
#' @param ... \code{\link{deep_mutational_scan}} objects to combine.
#' @return A \code{\link[tibble]{tibble}} whose rows contain positional data from all the provided scans.
#'
#' @export
rbind.deep_mutational_scan <- function(...) {
  df <- dplyr::bind_rows(lapply(list(...), as_tibble, full = TRUE))
  return(validate_combined_dms(df))
}

# TODO - check if scores are distributed properly or such?
#' Check a data frame is a combined mutational scan dataset
#'
#' Validate a data frame to determine if it sufficiently resembles a combined deep mutational scan dataset, as produced
#' by \code{\link{rbind.deep_mutational_scan}}. A error is raised if the data is incorrect, otherwise it is returned
#' as a \code{\link[tibble]{tibble}}.
#'
#' The following tests are performed:
#' \itemize{
#'   \item x is a \code{\link{data.frame}}
#'   \item There are no NA ER scores
#'   \item The required columns are present
#'   \item The expected additional columns are present for annotated data
#' }
#'
#' @param x Data frame to check.
#' @param annotated logical indicating the data be annotated with PCA, UMAP and Cluster information.
#' @return x converted to a \code{\link[tibble]{tibble}}
#' @export
validate_combined_dms <- function(x, annotated = FALSE) {
  if (!is.data.frame(x)) {
    stop("Combined mutational scan datasets must inherit from data.frame")
  }

  if (any(is.na(x[amino_acids]))) {
    stop("NA ER scores present\nUse impute_dms() to remove these before combining data")
  }

  cols <- c("study", "gene", "position", "wt", amino_acids)
  if (!all(cols %in% names(x))) {
    stop("Data frame does not contain required columns. Must contain: ", stringr::str_c(cols, collapse = ", "))
  }

  if (annotated) {
    cols <- c("cluster", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6",
              "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14",
              "PC15", "PC16", "PC17", "PC18", "PC19", "PC20", "umap1", "umap2",
              "base_cluster", "permissive", "small_margin", "high_distance",
              "dist1", "dist2", "dist3", "dist4", "dist5", "dist6", "dist7",
              "dist8", "notes")
    if (!all(cols %in% names(x))) {
      stop("Data frame does not contain the required columns. Annotated data must contain: ",
           stringr::str_c(cols, collapse = ", "))
    }
  }

  return(tibble::as_tibble(x))
}
