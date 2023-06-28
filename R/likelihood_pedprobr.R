#' Pedigree likelihood by pedprobr
#'
#' This uses [pedprobr::likelihood()] to compute the pedigree likelihood.
#'
#' @param x A `ped` object with at least one attached marker.
#' @param verbose A logical
#' @param ... Further arguments passed on to `pedprobr::likelihood()`
#'
#' @return A list with 3 entries:
#'
#'   * `program` : "pedprobr"
#'   * `likelihood` : the likelihood as computed by [pedprobr::likelihood()]
#'   * `time` : timing in seconds
#'
#' @export
likelihood_pedprobr = function(x, verbose = TRUE, ...) {
  if(verbose) cat("Program `pedprobr`...")
  if(!requireNamespace("pedprobr", quietly = TRUE)) {
    if(verbose) cat("skipped. Package not installed\n")
    return()
  }

  st = Sys.time()
  res = pedprobr::likelihood(x, markers = 1, verbose = FALSE, ...)
  time = format(round(Sys.time() - st, 2))
  if(verbose)
    cat(sprintf("finished in %s\n", time))

  list(program="pedprobr", likelihood=res, time=time)
}

