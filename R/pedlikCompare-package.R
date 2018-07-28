#' pedlikCompare: Compare pedigree likelihoods calculated by different programs
#'
#' Several R packages compute pedigree likelihoods, including pedtools,
#' paramlink, Familias, ElstonStewart. Outside of R, a widely used program is
#' MERLIN. This purpose of this package is to facilitate comparisons of these
#' programs, both in terms of numeric accuracy and timings. The pedtools package
#' is used as starting point for creatingpedigrees and markers, which are then
#' converted to other formats as needed.
#'
#' @docType package
#' @import pedtools
#' @import tibble
#'
#' @name pedlikCompare
NULL
