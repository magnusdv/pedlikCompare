#' Compare likelihood computations
#'
#' This is the main function of the package. For a given pedigree and a single
#' marker object, it computes the pedigree likelihood using different programs
#' and outputs a table containing the numerical results and timing. By default
#' the R packages being compared are `pedprobr`, `paramlink`, `Familias` and
#' `ElstonStewart`. In addition the external software `MERLIN` may be included
#' (via the wrapper pedprobr::merlin) if it is installed on the users computer.
#'
#' @param x A `pedtools::ped` object
#' @param marker Either a `pedtools::marker` object, or the name (or index) of
#'   an attached marker
#' @param verbose A logical
#' @param programs A character indicating which programs should be included. One
#'   of more of the terms "pedprobr", "paramlink", "merlin", "Familias", "ES".
#'   By default all are included.
#'
#' @references For MERLIN, see <http://csg.sph.umich.edu/abecasis/Merlin/>.

#' @export
compare = function(x, marker=1, verbose=TRUE,
                   programs=c("pedprobr", "paramlink", "Familias", "ES",
                              "merlin")) {
  if(is.marker(marker)) {
    x = setMarkers(x, list(marker))
    marker = 1
  }
  x = selectMarkers(x, marker)

  #if(is_Xmarker(x$markerdata[[1]])) stop("Sorry, only autosomal markers for now.")

  RESULT = tibble(program=character(), likelihood=numeric(), time=numeric())

  for(prog in programs) {
    FUN = get(sprintf("likelihood_%s", prog))
    res = FUN(x, verbose=verbose)
    if(length(res))
      RESULT = add_case(RESULT, !!!res)
  }

  RESULT = add_column(RESULT, lnlik = log(RESULT$likelihood), .after = "likelihood")
  RESULT
}
