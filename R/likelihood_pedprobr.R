#' Pedigree likelihood by pedprobr
#'
#' This uses [pedprobr::likelihood()] to compute the pedigree likelihood.
#'
#' @param x A `ped` object with at least one attached marker.
#' @param unit Unit for reporting runtimes, e.g. "auto" (default) or "secs".
#' @param verbose A logical.
#' @param ... Further arguments passed on to `pedprobr::likelihood()`.
#'
#' @return A list with 3 entries:
#'
#' * `program`: "pedprobr"
#' * `likelihood`: the likelihood computed by [pedprobr::likelihood()]
#' * `time`: runtime
#'
#' @export
likelihood_pedprobr = function(x, unit = "auto", verbose = TRUE, debug = FALSE, ...) {
  if(verbose) cat("Program `pedprobr`...")
  if(!requireNamespace("pedprobr", quietly = TRUE)) {
    if(verbose) cat("skipped. Package not installed\n")
    return()
  }

  st = Sys.time()
  res = pedprobr::likelihood(x, markers = 1, verbose = debug, ...)
  time = difftime(Sys.time(), st, units = unit)

  if(verbose)
    cat(sprintf("finished in %s\n", format(round(time, 2))))

  list(program="pedprobr", likelihood=res, time=time)
}

