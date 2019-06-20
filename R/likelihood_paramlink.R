#' Pedigree likelihood by paramlink
#'
#' Converts the `ped` object to a [paramlink::linkdat] object, and uses
#' [paramlink::likelihood()] to Compute the likelihood of the first marker.
#'
#' @param x A `ped` object
#' @param verbose A logical
#' @param ... Further arguments passed on to `paramlink::likelihood()`
#'
#' @return A list with 3 entries:
#'
#'   * `program` : "paramlink"
#'   * `likelihood` : the likelihood as computed by [paramlink::likelihood()]
#'   * `time` : timing in seconds
#'
#' @export
likelihood_paramlink = function(x, verbose=T, ...) {
  if(verbose) cat("Program `paramlink`...")
  if(!requireNamespace("paramlink", quietly = TRUE)) {
    if(verbose) cat("skipped. Package not installed\n")
    return()
  }

  y = ped2linkdat(x, verbose=F)

  st = Sys.time()
  res = paramlink::likelihood(y, locus1=1, ...)

  time = as.numeric(Sys.time() - st)
  if(verbose) cat(sprintf("finished in %.2f seconds\n", time))

  list(program="paramlink", likelihood=res, time=time)
}


ped2linkdat = function(x, verbose=F) {
  if (!requireNamespace("paramlink", quietly = TRUE))
    stop2("Package 'paramlink' is not installed")

  famid = x$FAMID
  if(famid == "") famid = 1

  mlist = x$markerdata

  x$markerdata = NULL
  p = cbind(famid, as.matrix(x), 1)
  colnames(p) = c("FAMID", "ID", "FID", "MID", "SEX", "AFF")

  y = paramlink::linkdat(p, verbose=verbose)

  if(!is.null(mlist)) {
    mlist = lapply(mlist, function(m) {
      attributes(m) =
        list(dim = dim(m),
             name = name(m),
             chrom = if(isXmarker(m)) 23 else chrom(m),
             pos = posCm(m),
             nalleles = nAlleles(m),
             alleles = alleles(m),
             afreq = as.vector(afreq(m)),
             missing = 0,
             mutmat = mutmod(m),
             class = "marker")
      m
    })
    y = paramlink::setMarkers(y, mlist)
  }
  y
}

