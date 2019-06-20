#' Pedigree likelihood by MERLIN
#'
#' This uses the MERLIN wrapper [pedprobr::likelihoodMerlin()] to Compute
#' the likelihood of the first marker. Skipped if the marker models mutations
#' (since MERLIN does not handle this).
#'
#' @param x A `ped` object with at least one attached marker.
#' @param verbose A logical
#'
#' @return A list with 3 entries:
#'
#'   * `program` : "merlin"
#'   * `likelihood` : the likelihood as computed by [pedprobr::likelihoodMerlin()]
#'   * `time` : timing in seconds
#'
#' @export
likelihood_merlin = function(x, verbose=T) {
  if(verbose) cat("Program `merlin`...")
  if(!requireNamespace("pedprobr", quietly = TRUE)) {
    if(verbose) cat("skipped. Package not installed\n")
    return()
  }
  if(allowsMutations(x$markerdata[[1]])) {
    if(verbose) cat("skipped. Marker allows mutations\n")
    return()
  }
  st = Sys.time()
  res = pedprobr::likelihoodMerlin(x, markers=1)

  time = as.numeric(Sys.time()-st)
  if(verbose) cat(sprintf("finished in %.2f seconds\n", time))

  list(program="merlin", likelihood=res, time=time)
}


