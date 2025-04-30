#' Compare likelihood computations
#'
#' This is the main function of the package. For a given pedigree and a single
#' marker object, it computes the pedigree likelihood using different programs
#' and outputs a table containing the numerical results and timing. By default
#' the R packages being compared are `pedprobr`, `Familias` and `ElstonStewart`.
#' In addition the external software `MERLIN` may be included (via the wrapper
#' pedprobr::merlin) if it is installed on the users computer.
#'
#' @param x A `pedtools::ped` object.
#' @param marker Either a `pedtools::marker` object, or the name (or index) of
#'   an attached marker.
#' @param theta Theta correction, by default 0. Of the supported packages, only
#'   "pedprobr" and "Familias" support theta correction.
#' @param unit Unit for reporting runtimes, e.g. "auto" (default) or "secs".
#' @param verbose A logical.
#' @param programs A character containing some of all of the words "pedprobr",
#'   "merlin", "Familias", "ES", indicating which programs should be included in
#'   the comparison. By default all except "ES" are included.
#' @param ... Further arguments passed on to [pedprobr::likelihood()].
#'
#' @references For MERLIN, see <https://csg.sph.umich.edu/abecasis/merlin/>.
#'
#' @examples
#' x = nuclearPed(2) |> addMarker(`3` = "1/2", `4` = "1/2")
#'
#' compare(x)
#'
#' compare(x, theta = 0.1)
#'
#' @importFrom cli symbol
#' @importFrom crayon green red white
#' @export
compare = function(x, marker = 1, theta = 0, unit = "auto", verbose = TRUE,
                   programs = c("pedprobr", "Familias", "merlin"), ...) {
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

  if(theta > 0) {
    rem = setdiff(programs, c("pedprobr", "Familias"))
    message("Ignoring programs without theta correction: ", toString(rem))
    programs = setdiff(programs, rem)
  }

  RESULT = tibble(program = character(), likelihood = numeric(), time = as.difftime(character()))

  for(prog in programs) {

    res = tryCatch(
      switch(prog,
             pedprobr = likelihood_pedprobr(x, theta = theta, unit = unit, verbose = verbose, ...),
             Familias = likelihood_Familias(x, kinship = theta, unit = unit, verbose = verbose),
             merlin = likelihood_merlin(x, unit = unit, verbose = verbose),
             ES = likelihood_ES(x, unit = unit, verbose = verbose)),
      error = function(e) e)

    if(inherits(res, "error")) {
      if(verbose) message(toString(conditionMessage(res)))
      next
    }

    if(length(res))
      RESULT = add_case(RESULT, !!!res)
  }

  RESULT = add_column(RESULT, lnlik = log(RESULT$likelihood), .after = "likelihood")

  if(verbose) {
    if(all_agree(RESULT))
      mess = green(cli::symbol$tick, "ALL PROGRAMS AGREE")
    else
      mess = red(cli::symbol$cross, "ANSWERS DISAGREE")
    cat("\n", mess, "\n\n")
  }

  RESULT$time = round(RESULT$time, 3)
  RESULT
}
