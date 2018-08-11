#' Pedigree likelihood by Familias
#'
#' Converts from `ped` to Familias objects (pedigree, loci, datamatrix) and uses
#' [Familias::FamiliasPosterior()] to Compute the pedigree likelihood.
#'
#' @param x A `ped` object
#' @param verbose A logical
#'
#' @return A list with 3 entries:
#'
#'   * `program` : "Familias"
#'   * `likelihood` : the likelihood as computed by [Familias::FamiliasPosterior()]
#'   * `time` : timing in seconds
#'
#' @export
likelihood_Familias = function(x, verbose=T) {
  if(verbose) cat("Program `Familias`...")
  if(!requireNamespace("Familias", quietly = TRUE)) {
    if(verbose) cat("skipped. Package not installed\n")
    return()
  }
  if(is_Xmarker(x$markerdata[[1]])) {
    if(verbose) cat("skipped. X-linked markers are not implemented for Familias\n")
    return()
  }

  FamiliasData = ped2Familias(x)

  st = Sys.time()
  res = Familias::FamiliasPosterior(FamiliasData$pedigree,
                                    FamiliasData$loci,
                                    FamiliasData$datamatrix)
  res = res$likelihoods

  time=as.numeric(Sys.time() - st)
  if(verbose) cat(sprintf("finished in %.2f seconds\n", time))

  list(program="Familias", likelihood=res, time=time)
}



#############################
### Familias conversions
#############################


ped2Familias = function(x) {
  pedigree = ped2FamiliasPedigree(x)
  loci = ped2FamiliasLoci(x)
  datamatrix = ped2FamiliasDatamatrix(x)

  list(pedigree=pedigree, loci=loci, datamatrix=datamatrix)
}

ped2FamiliasPedigree = function(x) {
  id = x$LABELS
  dadid = momid = rep(NA, pedsize(x))
  dadid[x$FID > 0] = id[x$FID]
  momid[x$MID > 0] = id[x$MID]
  sex = ifelse(x$SEX == 1, "male", "female")
  Familias::FamiliasPedigree(id, dadid, momid, sex)
}

ped2FamiliasDatamatrix = function(x) {
  if(!hasMarkers(x))
    return(NULL)

  ids = x$LABELS
  mlist = x$markerdata

  # TODO: lag alleleMatrix() i pedtools
  allelematrix = pedtools:::.prettyMarkers(mlist, missing=NA)

  as.data.frame(allelematrix, row.names=ids, stringsAsFactors=F)
}

ped2FamiliasLoci = function(x) {
  if(!hasMarkers(x))
    return(NULL)

  ## Replace NA marker names with dummy names
  #mnames = name(x, 1:nMarkers(x))
  #if(anyNA(mnames)) {
  #  idx = which(is.na(mnames))
  #  # name(x, idx) = paste0("m", idx) 
  #}

  lapply(x$markerdata, marker2FamiliasLocus)
}

marker2FamiliasLocus = function(m) {
  als = alleles(m)
  afr = afreq(m)
  mname = name(m)
  mutmat = attr(m, "mutmat")

  if(is.null(mutmat))
    Familias::FamiliasLocus(allelenames = als, frequencies = afr, name = mname)
  else
    Familias::FamiliasLocus(allelenames = als, frequencies = afr, name = mname,
                            MutationModel = "Custom",
                            maleMutationMatrix = mutmat$male,
                            femaleMutationMatrix = mutmat$female)
}
