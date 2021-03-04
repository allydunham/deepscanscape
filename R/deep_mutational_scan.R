# S3 class / attributes thing to store a dms dataset

# TODO - Store attributes to specify state or use list?
# TODO - data frame or tibble?
#' Deep mutational scan data
#'
#' Store deep mutational scanning data in a standardised format, alongside
#' metadata describing the state of the data. The primary data is contained in
#' a data frame
#'
#' @param x A data frame to parse
#' @param format Original data format (see description)
#' @return A deep mutational scan S3 object
#' @export
deep_mutational_scan <- function(x){

  class(x) <- c('deep_mutational_scan', class(x))
  return(x)
}

#' Print method for Deep Mutational Scans
#'
#' @export
print.deep_mutational_scan <- function(x, ...){

}

#' Summary method for Deep Mutational Scans
#'
#' @export
summary.deep_mutational_scan <- function(object, ...){

}

#' Autoplot method for Deep Mutational Scans
#'
#' @export
#' @importFrom ggplot2 autoplot
autoplot.deep_mutational_scan <- function(x, ...){

}

#' Plot method for Deep Mutational Scans
#'
#' @export
#' @importFrom graphics plot
plot.deep_mutational_scan <- function(x, ...){
  print(autoplot(x, ...))
}
