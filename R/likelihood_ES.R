#' Pedigree likelihood by ElstonStewart
#'
#' Converts the pedigree and its first marker to a [ElstonStewart::es.pedigree]
#' object and uses [ElstonStewart::Elston()] to Compute the pedigree likelihood.
#' This is only implemented for the case when the marker has exactly 2 alleles
#' and does not have a nontrivial mutation model.
#'
#' @param x A `ped` object with at least one attached marker.
#' @param verbose A logical
#'
#' @return A list with 3 entries:
#'
#'   * `program` : "ElstonStewart"
#'   * `likelihood` : the likelihood as computed by [ElstonStewart::Elston()]
#'   * `time` : timing in seconds
#'
#' @export
likelihood_ES = function(x, verbose=T) {
  if(verbose) cat("Program `ElstonStewart`...")
  if(!requireNamespace("ElstonStewart", quietly = TRUE)) {
    if(verbose) cat("skipped. Package not installed\n")
    return()
  }
  m = getMarkers(x, 1)[[1]]
  if(nAlleles(m) != 2) {
    if(verbose) cat("skipped. Non-diallelic markers are not implemented for the ElstonStewart package\n")
    return()
  }
  if(allowsMutations(m)) {
    if(verbose) cat("skipped. Marker allows mutations\n")
    return()
  }
  if(is_Xmarker(m)) {
    if(verbose) cat("skipped. X-linked markers are not implemented for the ElstonStewart package\n")
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
  es_ped = ElstonStewart::es.pedigree(id=x$ID, father=x$FID, mother=x$MID, sex=x$SEX,
                                      pheno=rep(0, N), geno=geno, famid="1")
  modele.di = ElstonStewart::modele.di
  st = Sys.time()
  res = ElstonStewart::Elston(es_ped, modele.di, list(p=p))$result

  time = as.numeric(Sys.time() - st)
  if(verbose) cat(sprintf("finished in %.2f seconds\n", time))

  list(program="ElstonStewart", likelihood=res, time=time)
}
