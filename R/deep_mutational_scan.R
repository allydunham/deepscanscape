# S3 class to store a dms dataset
# TODO - Test these funcs

#' Internal deep_mutational_scan constructor
#'
#' Create a new raw deep_mutational_scan object. This is expected to have no annotations and no imputation but have
#' had scores transformed and normalised onto the standard scale.
#'
#' @param df A data frame in the required format (position, wt, A-Y, additional data...).
#' @param study The source study
#' @param gene The studied gene
new_deep_mutational_scan <- function(df, study, gene) {
  out <- list(data = tibble::as_tibble(df), study = study, gene = gene,
              annotated = FALSE, impute_mask = NA, cluster = NA)
  class(out) <- c("deep_mutational_scan")
  return(out)
}

# TODO Check position / wt rows are distinct, this is assumed downstream
#' Validate deep_mutational_scan objects
#'
#' @param x \link{deep_mutational_scan} object
validate_deep_mutational_scan <- function(x) {
  if (!"deep_mutational_scan" %in% class(x)) {
    stop("deep_mutational_scan not listed in class()")
  }

  if (!all(c("data", "study", "gene", "annotated", "impute_mask", "cluster") %in% names(x))) {
    stop("Missing field(s). The object must contain 'data', 'study', 'gene', 'annotated' & 'impute_mask'")
  }

  if (!is.logical(x$annotated)) {
    stop("Invalid 'annotated' value (", x$annotated, "). Must be TRUE/FALSE")
  }

  if (!(is.matrix(x$impute_mask) & is.logical(x$impute_mask) | is.na(x$impute_mask))) {
    stop("Invalid 'impute_mask'. Must be a logical matrix or NA")
  }

  headers <- c("position", "wt", amino_acids)[which(!c("position", "wt", amino_acids) %in% names(x$data))]
  if (length(headers) > 0) {
    stop("Missing data columns: ", paste(headers, sep = ", "))
  }

  # TODO - Check annotations present
  # TODO - Check scores fit right profile? (Maybe too fragile)
  return(x)
}

#' Deep mutational scan data
#'
#' Store deep mutational scanning data in a standardised format, alongside
#' metadata describing the state of the data. The primary data is contained in
#' a tibble
#'
#' @param df A data frame to parse
#' @param scheme Original data scheme (see \code{\link{parse_dms_data}})
#' @param trans Function to transform scores onto the standard scale. Accepts a string corresponding to known transforms
#' or a custom function (\code{\link{transform_dms}}).
#' @param na_value Value to set missense NA scores to (see \code{\link{impute_dms}})
#' @param annotate Annotate the dataset with mutational landscape data (PCA, UMAP and amino acid subtypes)
#' @param study Source of the deep mutational scan
#' @param gene Gene scanned
#' @return A deep mutational scan S3 object
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
  df <- dplyr::select(df, .data$position, .data$wt, amino_acids, dplyr::everything())

  # Construct object
  out <- validate_deep_mutational_scan(new_deep_mutational_scan(df = df, study = study, gene = gene))

  # Impute
  if (!(is.na(na_value) | is.null(na_value))) {
    out <- impute_dms(out)
  }

  # Annotate
  if (annotate) {
    out <- annotate_dms(out)
  }

  return(out)
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
#' @param x \link{deep_mutational_scan} object
#' @param i,j,... Indices to access
#' @param drop Coerce result to lowest possible dimension
#' @param value Value to set
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
#' @param x,object \link{deep_mutational_scan} object
#' @param indent.str Indent string for \link[utils]{str}
#' @param nest.lev Nest level for \link[utils]{str}
#' @param ... Additional arguments
#'
#' @name dms_s3
NULL
#> NULL

# TODO - Pretty print when extra columns added
#' @describeIn dms_s3 S3 print method
#' @export
format.deep_mutational_scan <- function(x, ...) {
  # Select columns to show
  # TODO add more columns when annotated?
  # TODO show additional added columns

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
  tbl <- dplyr::select(tbl, dplyr::one_of("cluster"), .data$position, .data$wt, amino_acids, dplyr::everything())
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

#' Convert deep_mutational_scans to data frames
#'
#' @param x \link{deep_mutational_scan} object
#' @param full Include columns containing metadata (study, gene)
#' @param ... Ignored
#'
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

#' Combine deep mutational scan data
#'
#' @param ... \link{deep_mutational_scan} objects to combine
#'
#' @export
rbind.deep_mutational_scan <- function(...) {
  dplyr::bind_rows(lapply(list(...), as_tibble, full = TRUE))
}

# TODO - Need output object for this?
#' Summarise Deep Mutational Scans
#'
#' @param object A \link{deep_mutational_scan} object.
#' @param ... Additional arguments
#'
#' @export
summary.deep_mutational_scan <- function(object, ...) {
  stop("Not implemented yet")
}

#' Summary plot for Deep Mutational Scans
#'
#' @param object,x \link{deep_mutational_scan} to produce a
#' summary plot for.
#' @param ... Additional parameters
#'
#' @export
#' @name dms_plot
#' @importFrom ggplot2 autoplot
autoplot.deep_mutational_scan <- function(object, ...){ # nolint
  stop("Not implemented yet")
}

#' @describeIn dms_plot S3 plot method
#' @export
#' @importFrom graphics plot
plot.deep_mutational_scan <- function(x, ...) {
  print(autoplot(x, ...))
}
