#' Pedigree likelihood by Familias
#'
#' Converts from `ped` to Familias objects (pedigree, loci, datamatrix) and uses
#' [Familias::FamiliasPosterior()] to Compute the pedigree likelihood.
#'
#' @param x A `ped` object with at least one attached marker.
#' @param unit Unit for reporting runtimes, e.g. "auto" (default) or "secs".
#' @param verbose A logical.
#' @param ... Further arguments passed on to `FamiliasPosterior()`
#'
#'
#' @return A list with 3 entries:
#'
#' * `program`: "Familias"
#' * `likelihood`: the likelihood computed by [Familias::FamiliasPosterior()]
#' * `time`: runtime
#'
#' @export
likelihood_Familias = function(x, unit = "auto", verbose = TRUE, ...) {
  if(verbose) cat("Program `Familias`...")
  if(!requireNamespace("Familias", quietly = TRUE)) {
    if(verbose) cat("skipped. Package not installed\n")
    return()
  }
  if(isXmarker(x, 1)) {
    if(verbose) cat("skipped. X-linked markers are not implemented for Familias\n")
    return()
  }

  FamiliasData = ped2Familias(x)

  st = Sys.time()
  res = Familias::FamiliasPosterior(FamiliasData$pedigree,
                                    FamiliasData$loci,
                                    FamiliasData$datamatrix,
                                    ...)
  res = res$likelihoods
  time = difftime(Sys.time(), st, units = unit)

  if(verbose)
    cat(sprintf("finished in %s\n", format(round(time, 2))))


  list(program = "Familias", likelihood = res, time = time)
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
  id = labels(x)
  dadid = momid = rep(NA, pedsize(x))
  dadid[x$FIDX > 0] = id[x$FIDX]
  momid[x$MIDX > 0] = id[x$MIDX]
  sex = ifelse(x$SEX == 1, "male", "female")
  Familias::FamiliasPedigree(id, dadid, momid, sex)
}

ped2FamiliasDatamatrix = function(x) {
  getAlleles(x)
}

ped2FamiliasLoci = function(x) {
  if(!hasMarkers(x))
    return(NULL)

  lapply(x$MARKERS, marker2FamiliasLocus)
}

#' @importFrom tools toTitleCase
marker2FamiliasLocus = function(m) {
  als = alleles(m)
  afr = afreq(m)
  mname = name(m)
  mutmod = attr(m, "mutmod")

  args = list(allelenames = als, frequencies = afr, name = mname)

  if(!is.null(mutmod)) {
    params = pedmut::getParams(mutmod, format = 2)

    if(params[["model.F"]] %in% c("equal", "proportional", "stepwise")) {
      args$femaleMutationModel = params[["model.F"]] |> tools::toTitleCase()
      args$femaleMutationRate  = params[["rate.F"]]
      args$femaleMutationRate2  = params[["rate2.F"]]
      args$femaleMutationRange  = params[["range.F"]]
    }
    else {
      args$femaleMutationModel = "Custom"
      args$femaleMutationMatrix = mutmod$female
    }
    if(params[["model.M"]] %in% c("equal", "proportional", "stepwise")) {
      args$maleMutationModel = params[["model.M"]] |> tools::toTitleCase()
      args$maleMutationRate  = params[["rate.M"]]
      args$maleMutationRate2  = params[["rate2.M"]]
      args$maleMutationRange  = params[["range.M"]]
    }
    else {
      args$maleMutationModel = "Custom"
      args$maleMutationMatrix = mutmod$male
    }
  }

  do.call(Familias::FamiliasLocus, args)
}

