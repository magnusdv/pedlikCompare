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
#' @references For MERLIN, see <https://csg.sph.umich.edu/abecasis/merlin/>.
#'
#' @importFrom crayon bgGreen bgRed bold
#' @export
compare = function(x, marker=1, verbose=TRUE,
                   programs=c("pedprobr", "paramlink", "Familias", "ES",
                              "merlin")) {
  if(!is.ped(x))
    stop2("Input is not a `ped` object")
  if(is.marker(marker)) {
    x = setMarkers(x, list(marker))
    marker = 1
  }
  else {
    if(!hasMarkers(x)) stop2("The pedigree has no attached markers")
    midx = whichMarkers(x, marker)
    if(length(midx) == 0) stop2("Marker not found")
    if(length(midx) > 1) stop2("Multiple markers selected")
    x = selectMarkers(x, midx)
  }

  programs = match.arg(programs, several.ok = TRUE)

  RESULT = tibble(program=character(), likelihood=numeric(), time=character())
  for(prog in programs) {

    FUN = get(sprintf("likelihood_%s", prog))

    res = tryCatch(FUN(x, verbose=verbose), error=function(e) e)
    if(inherits(res, "error")) {
      if(verbose) message(toString(res))
      next
    }

    if(length(res))
      RESULT = add_case(RESULT, !!!res)
  }

  RESULT = add_column(RESULT, lnlik = log(RESULT$likelihood), .after = "likelihood")

  if(verbose) {
    check  = all_agree(RESULT)
    if(check)
      cat(crayon::bgGreen$bold("===> ALL PROGRAMS AGREE! <===\n"))
    else
      cat(crayon::bgRed$bold("===> ANSWERS ARE NOT THE SAME <===\n"))
  }

  RESULT
}
