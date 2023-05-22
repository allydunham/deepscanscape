# S3 class to store a dms dataset

#' Internal deep_mutational_scan constructor
#'
#' Create a new raw deep_mutational_scan object, containing data from one or more studies.
#' This is expected to have no annotations and no imputation but have had all scores transformed and
#' normalised onto the standard scale.
#'
#' @param df A data frame in the required format (name, position, wt, A-Y, additional data...).
#' @param meta A data frame detailing the metadata for each study. Each row should contain the name, study, gene,
#' source, description and positions entries, covering each unique name in df.
#' @return A deep_mutational_scan S3 object (see \code{\link{deep_mutational_scan}})
new_deep_mutational_scan <- function(df, meta) {
  out <- list(data = tibble::as_tibble(df), meta = meta, annotated = FALSE,
              imputed = FALSE, multi_study = nrow(meta) > 1)
  class(out) <- c("deep_mutational_scan")
  return(out)
}

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
#'   \item \code{imputed} is logical
#'   \item \code{annotated} is logical
#'   \item \code{multi_study} is logical
#'   \item \code{data} contains the correct columns and no incorrect ones
#'   \item \code{data} contains no duplicate rows
#'   \item \code{meta} list is correctly structured
#'   \item \code{meta} list contains values for each unique name
#'   \item The correct additional columns are present when the object is annotated or imputed
#'   \item No NA scores remain for imputed objects
#'   \item Objects listed as multi_study contain data from multiple studies
#' }
#'
#' @param x \code{\link{deep_mutational_scan}} object.
#' @return The input is returned unaltered, since the main purpose of the function is to raise an error if the
#'   object is erroneous.
validate_deep_mutational_scan <- function(x) { # nolint
  # Expected headers
  data_headers <- c("name", "position", "wt", amino_acids, "stop")
  impute_headers <- paste0("impute_", amino_acids)
  annotated_headers <- c(paste0("PC", seq_len(20)), "umap1", "umap2", "cluster", "base_cluster", "permissive",
                         "ambiguous", "high_distance", paste0("dist", seq_len(8)), "cluster_notes")
  meta_headers <- c("name", "study", "gene", "source", "positions")

  # Check fundamental structure
  if (!is.deep_mutational_scan(x)) {
    stop("Not a deep_mutational_scan (check class(x))")
  }

  if (!all(c("data", "meta", "annotated", "imputed") %in% names(x))) {
    stop("Missing field(s). The object must contain 'data', 'meta', 'annotated' & 'imputed'")
  }

  # Check field types
  if (!tibble::is_tibble(x$data)) {
    stop("Invalid 'data' object. Must be a tibble")
  }

  if (!tibble::is_tibble(x$meta)) {
    stop("Invalid 'meta' object. Must be a tibble")
  }

  if (!is.logical(x$imputed)) {
    stop("Invalid 'imputed' value (", x$imputed, "). Must be TRUE/FALSE")
  }

  if (!is.logical(x$annotated)) {
    stop("Invalid 'annotated' value (", x$annotated, "). Must be TRUE/FALSE")
  }

  if (!is.logical(x$multi_study)) {
    stop("Invalid 'multi_study' value (", x$multi_study, "). Must be TRUE/FALSE")
  }

  # Check data object
  missing_headers <- data_headers[which(!data_headers %in% names(x$data))]
  if (length(missing_headers) > 0) {
    stop("Missing data columns: ", paste(missing_headers, sep = ", "))
  }

  unexpected_headers <- names(x$data)[!names(x$data) %in% c(data_headers, impute_headers, annotated_headers)]
  if (length(unexpected_headers) > 0) {
    stop("Data contains unsupported columns:", paste(unexpected_headers, sep = ", "))
  }

  if (nrow(x$data) != nrow(dplyr::distinct(x$data[c("name", "position", "wt")]))) {
    stop("Data contains duplicated positions")
  }

  # Check meta object
  missing_headers <- meta_headers[which(!meta_headers %in% names(x$meta))]
  if (length(missing_headers) > 0) {
    stop("Missing meta columns: ", paste(missing_headers, sep = ", "))
  }

  if (!all(unique(x$data$name) %in% unique(x$meta$name))) {
    stop("Names in meta and data don't match")
  }

  # Check imputation
  if (x$imputed) {
    if (any(is.na(x$data[amino_acids]))) {
      stop("Scan is marked as imputed by has NA ER scores")
    }

    missing_headers <- impute_headers[which(!impute_headers %in% names(x$data))]
    if (length(missing_headers) > 0) {
      stop("Scan is marked as imputed by is missing data columns: ", paste(missing_headers, sep = ", "))
    }

    if (any(!as.matrix(x$data[impute_headers]) %in% c(0, 1, 2))) {
      stop("Unexpected values in impute_X columns. All should be 0, 1 or 2")
    }
  }

  # Check annotation
  if (x$annotated) {
    missing_headers <- annotated_headers[which(!annotated_headers %in% names(x$data))]
    if (length(missing_headers) > 0) {
      stop("Scan is marked as annotated by is missing data columns: ", paste(missing_headers, sep = ", "))
    }
  }

  return(x)
}

#' Deep mutational scan data
#'
#' Store deep mutational scanning data in a standardised format, alongside
#' metadata describing the the data. This function creates a deep_mutational_scan
#' S3 object. Several of these can be joined using \code{\link{bind_scans}} to create a
#' combined dataset.
#'
#' @param df A data frame to parse.
#' @param name Name of the deep mutational scan.
#' @param scheme Original data scheme (see \code{\link{parse_deep_scan}}).
#' @param trans Function to transform scores onto the standard scale. Accepts a string corresponding to known transforms
#' or a custom function (\code{\link{transform_er}}).
#' @param na_value How to set missense NA scores (see \code{\link{impute}}).
#' @param annotate Annotate the dataset with mutational landscape data (PCA, UMAP and amino acid subtypes).
#' @param study Study in which the scan was performed.
#' @param gene Gene scanned.
#' @param source Source of study. This is for reference only and is not used internally.
#' @param description Description of study, containing any miscellaneous details. This is for reference only and is
#' not used internally.
#' @param ... Additional arguments passed to \code{\link{parse_deep_scan}}.
#' @return A deep_mutational_scan S3 object, containing the following fields:
#' \itemize{
#'   \item data: wide format \code{\link[tibble]{tibble}} containing ER scores and other positional data
#'   \item meta: Tibble containing meta data about each study in the dataset
#'   \item imputed: logical indicating the dataset has been completed via imputation(see \code{\link{impute}}
#'   for details)
#'   \item annotated: logical indicating if the data has been annotated with PCs, UMAP coordinates and clusters
#'   \item multi_study: logical indicating the dataset contains data from multiple studies
#' }
#'
#' The \code{data} tibble contains the following fields (those marked * are not always present):
#' \itemize{
#'   \item name: Name of the scan this position is from.
#'   \item position: Position in the protein.
#'   \item wt: Wild type amino acid.
#'   \item A-Y: Fitness for each substitution.
#'   \item stop: Fitness for early stop substitutions
#'   \item impute_A-impute_Y: Whether the fitness score is imputed (see \code{\link{impute}}).
#'   \item PC1-PC20: Deep landscape principal component coordinates.
#'   \item umap1-umap2: Deep landscape UMAP coordinates.
#'   \item cluster: Assigned amino acid subtype.
#'   \item base_cluster: Nearest functional subtypes centroid in PC cosine space, which is the assigned subtype before
#'   corrections are made for permissive, outlier and ambiguous subtypes.
#'   \item permissive, ambiguous, high_distance: Whether the position fulfilled the criteria for the special subtypes.
#'   \item dist1-dist8: PC cosine distance to each cluster centroid for subtypes of that amino acid.
#'   \item cluster_notes: Notes about cluster assignment.
#' }
#'
#' The class is used like a list apart from [ accesses the main data tibble (see \link[=dms_extract]{details}).
#' The class is also associated with other common generics (see \link[=dms_s3]{basic functions},
#' \link[=dms_summary]{summary}, \link[=as_tibble]{data frames}, \link[=dms_plot]{plotting}. Multiple deep_mutational
#' scan objects can be combined using \link{bind_scans}).
#' @examples
#'
#' # Create a scan object
#' path <- system.file("extdata", "urn_mavedb_00000011_a_1_scores.csv",
#'                     package = "deepscanscape")
#' csv <- read.csv(path, skip = 4)
#' dms <- deep_mutational_scan(csv, name = "Hietpas Hsp90", scheme = "mave", trans = NULL,
#'                             na_value = "impute", annotate = FALSE, gene = "Hsp90",
#'                             study = "Hietpas et al. (2011)", source = "",
#'                             description = "Scan of 9 hsp90 positions")
#' @export
deep_mutational_scan <- function(df, name, scheme=NULL, trans=NULL, na_value="impute",
                                 annotate=FALSE, study=NA, gene=NA, source=NA, description=NA, ...) {
  # Parse scheme
  if (!is.null(scheme)) {
    df <- parse_deep_scan(df, scheme = scheme, ...)
  } else {
    df <- tibble::as_tibble(df)
    if (!all(c("position", "wt", "mut", "score") %in% names(df))) {
      stop("Incorrect input df, must include position, wt, mut and score columns")
    }
    df <- dplyr::select(df, .data$position, .data$wt, .data$mut, .data$score)
  }

  # Transform
  if (!is.null(trans)) {
    df$score <- transform_er(df$score, trans)
  }

  # Normalise
  df$score <- normalise_er(df$score)

  # Format tibble
  df <- df[df$mut %in% c(amino_acids, "stop"), ]
  df <- dplyr::arrange(df, .data$position, .data$mut)
  df <- tidyr::pivot_wider(df, names_from = "mut", values_from = "score")
  exp_cols <- structure(rep(NA_real_, 21), names = c(amino_acids, "stop"))
  df <- tibble::add_column(df, !!!exp_cols[setdiff(names(exp_cols), names(df))])
  df$name <- name
  df <- dplyr::select(df, .data$name, .data$position, .data$wt, dplyr::all_of(amino_acids), .data$stop,
                      dplyr::everything())

  # Construct metadata
  meta <- tibble::tibble(name = name, study = study, gene = gene, source = source,
                         description = description, positions = dplyr::n_distinct(df$position))

  # Construct object
  out <- validate_deep_mutational_scan(new_deep_mutational_scan(df = df, meta = meta))

  # Impute
  if (is.null(na_value)) {
    na_value <- NA
  }

  if (!is.na(na_value)) {
    out <- impute(out, na_value = na_value)
  } else {
    if (any(is.na(out$data[amino_acids]))) {
      warning("No imputation applied but NA values present. ",
              "Data may not be suitable for downstream analysis until NA values are removed.")
    }
  }

  # Annotate
  if (annotate) {
    out <- annotate(out)
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
#' dms <- deepscanscape::deep_scans$p53
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

#' @describeIn dms_s3 S3 format method
#' @export
format.deep_mutational_scan <- function(x, ...) {
  out <- "# Deep mutational scanning data"

  if (x$multi_study) {
    out <- c(out, stringr::str_c("# ", nrow(x$meta), " scans:"),
             stringr::str_c("# ", x$meta$study, " - ", x$meta$gene, " (", x$meta$positions,
                            " positions, name: ", x$meta$name, ")"))
  } else {
    out <- c(out,
             stringr::str_c("# Name: ", x$meta$name),
             stringr::str_c("# Study: ", x$meta$study),
             stringr::str_c("# Gene: ", x$meta$gene),
             stringr::str_c("# ", x$meta$positions, " positions"))
  }

  if (x$annotated) {
    out <- c(out, "# Annotated with clusters, pricipal components and UMAP coordinates")
  }

  out <- c(out, "# Positional data:")

  if (requireNamespace("crayon", quietly = TRUE)) {
    out <- crayon::make_style("darkgrey")(out)
  }



  cols <- c("position", "wt", amino_acids, "stop")
  if (x$annotated) {
    cols <- c("cluster", cols)
  }

  if (x$multi_study) {
    cols <- c("name", cols)
  }

  tbl <- dplyr::select(x$data, dplyr::all_of(cols), dplyr::everything())
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

  cat("deep_mutational_scan [", nrow(object$meta), " scans with ", sum(object$meta$positions), " total positions]",
      "\n", sep = "")

  utils::str(list(meta = object$meta, annotated = object$annotated, imputed = object$imputed,
                  multi_study = object$multi_study, data = as.list(object$data)),
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
#' @param ... Ignored.
#' @return A \code{\link[tibble]{tibble}}, \code{\link[base]{data.frame}} or \code{\link[base]{list}}
#' @importFrom tibble as_tibble
#' @name as_tibble
#' @export
as_tibble.deep_mutational_scan <- function(x, ...) { # nolint
  return(x$data)
}

#' @describeIn as_tibble S3 as.data.frame method
#' @export
as.data.frame.deep_mutational_scan <- function(x, ...) {
  return(as.data.frame(as_tibble(x, ...)))
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
  if (object$multi_study) {
    plot_landscape(object)
  } else {
    plot_er_heatmap(object)
  }
}

#' @describeIn dms_plot S3 plot method
#' @export
#' @importFrom graphics plot
plot.deep_mutational_scan <- function(x, ...) {
  print(autoplot(x, ...))
}

#' Combine deep mutational scan data
#'
#' Combine multiple \code{\link{deep_mutational_scan}} objects into a single \code{\link[tibble]{tibble}}, including
#' columns giving the study and gene of each position.
#'
#' @param ... \code{\link{deep_mutational_scan}} objects to combine or lists of such objects.
#' @param impute_missing Impute missing ER scores using \code{\link{impute}} where scans have not already been imputed.
#' @param annotate_missing Annotate unannotated scans using \code{\link{annotate}}.
#' @return A multi study \code{\link{deep_mutational_scan}}.
#'
#' @export
bind_scans <- function(..., impute_missing = FALSE, annotate_missing = FALSE) {
  scans <- rlang::flatten(list(...))

  if (!all(sapply(scans, is.deep_mutational_scan))) {
    stop("All objects must be deep_mutational_scan objects or lists of these")
  }

  imputed <- sapply(scans, `[[`, "imputed")
  annotated <- sapply(scans, `[[`, "annotated")

  imp <- all(imputed)
  ann <- all(annotated)

  if (impute_missing) {
    scans <- lapply(scans, function(x) if (x$imputed) x else impute(x))
    imp <- TRUE
  } else {
    if (any(imputed) & !imp) {
      stop("Mixture of imputed and unimputed scans. All must be imputed or all not imputed. ",
           "Use impute_missing = TRUE to impute ER values for unimputed scans.")
    }
  }

  if (annotate_missing) {
    scans <- lapply(scans, function(x) if (x$annotated) x else annotate(x))
    ann <- TRUE
  } else {
    if (any(imputed) & !imp) {
      stop("Mixture of annotated and unannotated scans. All must be annoteated or all not annotated. ",
           "Use annotate_missing = TRUE to annotate scans missing annotation.")
    }
  }

  scan_names <- sapply(scans, function(x) x$meta$name)

  if (anyDuplicated(scan_names)) {
    stop("Scans contain duplicate names. All scans must be uniquely identified via the 'name' column in data and meta.")
  }

  df <- dplyr::bind_rows(lapply(scans, `[[`, "data"))
  meta <- dplyr::bind_rows(lapply(scans, `[[`, "meta"))
  out <- new_deep_mutational_scan(df, meta)
  out$annotated <- ann
  out$imputed <- imp
  out$multi_study <- TRUE

  return(validate_deep_mutational_scan(out))
}
