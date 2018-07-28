#' Pedigree likelihood by pedprobr
#'
#' This uses [pedprobr::likelihood()] to Compute the pedigree likelihood.
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
likelihood_pedprobr = function(x, verbose=T, ...) {
  if(verbose) cat("Program `pedprobr`...")
  if(!requireNamespace("pedprobr", quietly = TRUE)) {
    if(verbose) cat("skipped. Package not installed\n")
    return()
  }

  st = Sys.time()
  res = pedprobr::likelihood(x, marker1=1, verbose=F, ...)
  time = as.numeric(Sys.time()-st)

  if(verbose) cat(sprintf("finished in %.2f seconds\n", time))
  list(program="pedprobr", likelihood=res, time=time)
}

