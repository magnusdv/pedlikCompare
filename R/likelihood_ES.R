#' Pedigree likelihood by ElstonStewart
#'
#' Converts the pedigree and its first marker to a
#' `ElstonStewart::es.pedigree()` object and uses `ElstonStewart::Elston()` to
#' compute the pedigree likelihood. This is only implemented for the case when
#' the marker has exactly 2 alleles and does not allow mutation models.
#'
#' @param x A `ped` object with at least one attached marker.
#' @param unit Unit for reporting runtimes, e.g. "auto" (default) or "secs".
#' @param verbose A logical.
#'
#' @return A list with 3 entries:
#'
#'   * `program` : "ElstonStewart"
#'   * `likelihood` : the likelihood as computed by `ElstonStewart::Elston()`
#'   * `time` : timing in seconds
#'
#' @export
likelihood_ES = function(x, unit = "auto", verbose=T) {
  if(verbose) cat("Program `ElstonStewart`...")
  if(!requireNamespace("ElstonStewart", quietly = TRUE)) {
    if(verbose) cat("skipped. Package not installed\n")
    return()
  }
  m = getMarkers(x, 1)[[1]]
  if(nAlleles(m) != 2) {
    if(verbose) cat("skipped (non-diallelic markers are not implemented)\n")
    return()
  }
  if(allowsMutations(m)) {
    if(verbose) cat("skipped (mutations are not implemented)\n")
    return()
  }
  if(isXmarker(m)) {
    if(verbose) cat("skipped (X-linked markers are not implemented)\n")
    return()
  }
  N = pedsize(x)
  ALTcount = rowSums(m==2)
  MISScount = rowSums(m==0)
  geno = rep(list(0:2), N)
  geno[ALTcount==2] = 2
  geno[ALTcount==1 & MISScount==0] = 1
  geno[ALTcount==0 & MISScount==0] = 0
  geno[ALTcount==1 & MISScount==1] = 1:2 # TODO Warning message: number of items to replace is not a multiple of replacement length
  geno[ALTcount==0 & MISScount==1] = 0:1

  p = afreq(m)[1]
  es_ped = ElstonStewart::es.pedigree(id=1:pedsize(x), father=x$FIDX, mother=x$MIDX, sex=x$SEX,
                                      pheno=rep(0, N), geno=geno, famid="1")
  modele.di = ElstonStewart::modele.di
  st = Sys.time()
  res = ElstonStewart::Elston(es_ped, modele.di, list(p=p))$result
  res = unname(res)   # remove annoying name

  time = difftime(Sys.time(), st, units = unit)

  if(verbose)
    cat(sprintf("finished in %s\n", format(round(time, 2))))

  list(program="ElstonStewart", likelihood=res, time=time)
}
