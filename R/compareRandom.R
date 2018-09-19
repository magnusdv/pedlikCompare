
#' @importFrom stats runif rpois
randomTestCase = function(pedname = NULL, subsettype = NULL, ids = NULL, swapsex = NULL,
                          Xchrom = NULL, nall = NULL, alsChar = NULL, alsPerm = NULL,
                          muttype = NULL, relab = NULL, reorder = NULL) {

  if(!requireNamespace("forrel", quietly = TRUE))
    stop2("Package `forrel` is not installed")
  if(!requireNamespace("pedmut", quietly = TRUE))
    stop2("Package `pedmut` is not installed")

  # Pedigrees to choose from
  PEDS = list(SING = singleton(1),
              NUC = nuclearPed(1),
              LIN = linearPed(2),
              HS = halfSibPed(),
              ANC = relabel(addParents(linearPed(2), 4, verbose = F), 1:7),
              FMS = fullSibMating(1),
              RAND = "random")

  # Specially interesting subsets of genotyped individuals.
  SUBSETS = list(SING = list(),
                 NUC = list(1,2,3,2:3),
                 LIN = list(1, 3, 5),
                 HS = list(4:5),
                 ANC = list(1,3,4,6,7, c(1,2), c(1,7), c(4,7), c(1,4,7)),
                 FMS = list(6, 5:6, c(1,6)),
                 RAND = list())

  # Pedigree
  if(is.null(pedname))
    pedname = sample(names(PEDS), size = 1)

  if(pedname == "RAND") {
    founders = rpois(1, 2) + 1
    g = sample(founders + 2:3, size= 1)
    ped = randomPed(g, founders)
    if(is.pedList(ped)) ped = ped[[1]]
  }
  else
    ped = PEDS[[pedname]]

  # Subset of genotyped individuals
  if(is.null(ids)) {
    randomsub = labels(ped)[sample(c(T,F), pedsize(ped), replace = T)]
    subsetList = c(Fixed = SUBSETS[[pedname]],
                   list(None = character(0),
                        Leaves = leaves(ped),
                        All = labels(ped),
                        Random = randomsub))

    if(is.null(subsettype))
      subsettype = sample(names(subsetList), size = 1)

    ids = subsetList[[subsettype]]
  }

  # permute genders?
  if(is.null(swapsex))
    swapsex = sample(c(T,F), 1)
  if(swapsex) {
    swap_ids = labels(ped)[sample(c(T,F), pedsize(ped), replace = T)]
    ped = swapSex(ped, swap_ids, verbose = F)
  }

  ######################
  ### Marker attributes


  # Number of alleles
  if(is.null(nall))
    nall = rpois(1, lambda=2) + 1 # ensure > 0

  # Alleles
  als = 1:nall

  # Permute alleles?
  if(is.null(alsChar))
    alsChar = sample(c(T,F), 1)
  if(alsChar) als = letters[als]

  # Permute alleles?
  if(is.null(alsPerm))
    alsPerm = sample(c(T,F), 1)
  if(alsPerm) als = sample(als)

  # Frequencies
  freqs = runif(nall, 0.01, 0.99)
  freqs = round(freqs/sum(freqs), 2)
  if(nall > 1) freqs[nall] = 1 - sum(freqs[-nall])

  # Mutation model
  # 0 = no model; 1 = samle male/female; 2 = different male/female
  if(is.null(muttype))
    muttype = sample(0:2, 1)

  models = sample(c("eq", "prop", "ran", "triv"), size = muttype)

  if(muttype > 0) {
    rate = round(runif(muttype, 0.01, 0.04), 2)
    if(muttype == 2) {
      models = structure(as.list(models), names = c("male", "female"))
      rate = structure(as.list(rate), names = c("male", "female"))
    }
    mutmod = pedmut::mutationModel(models, alleles = als, afreq=freqs, rate = rate)
  }
  else {
    mutmod = NULL
    rate = NULL
  }

  # X chromosome?
  if(is.null(Xchrom))
    Xchrom = (muttype == 0) && sample(c(T,F), 1)

  #####################
  #### Simulate 1 marker
  x = forrel::simpleSim(ped, N = 1, alleles = als, afreq = freqs, ids = ids,
                        mutmod = mutmod, Xchrom = Xchrom, verbose = F)

  m = x$markerdata[[1]]
  geno = paste(format(m)[ids], collapse = ", ")

  # Relabel?
  if(is.null(relab))
    relab = sample(c(T,F), 1)
  if(relab)
    x = relabel(x, sample(letters[1:pedsize(x)]))

  # Reorder?
  if(is.null(reorder))
    reorder = sample(c(T,F), 1)
  if(reorder)
    x = reorderPed(x, neworder = sample(1:pedsize(x)))

  list(ped = x, originalPed = ped, pedname = pedname, swapsex = swapsex,
       ids = paste(ids, collapse=","), Xchrom = Xchrom, geno = geno,
       nall = nall, alleles = paste(als, collapse=","), freqs = freqs,
       alsChar = alsChar, alsPerm = alsPerm,
       muttype = muttype, models = paste(unlist(models), collapse=","),
       mutrate = paste(unlist(rate), collapse=","),
       relabel = relab, reorder = reorder)
}

#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
compareRandom = function(n = 1, programs = NA, verbose = F, store_bad = FALSE, ...) {
  if(store_bad)
    BAD = list()

  # Progress bar
  useProgBar = n > 9 && !verbose
  if(useProgBar) pb = txtProgressBar(min = 0, max = n, style = 3)

  for(i in 1:n) {
    case = randomTestCase(...)

    if(is.na(programs))
      use_programs = c("pedprobr", ifelse(case$Xchrom, "merlin", "Familias"))
    else
      use_programs = programs

    if(verbose) {
      cat(sprintf("%d. %s, %d alleles. Mutmodel: %s.", i, case$pedname, case$nall, case$models))
      st = proc.time()
    }

    # Run likelihood calculations
    setTimeLimit(elapsed = 10, transient = TRUE)
    results = tryCatch({
      compare(case$ped, 1, programs = use_programs, verbose = F)
      }, error = function(e) NULL)
    setTimeLimit(elapsed = Inf, transient = TRUE)

    if(verbose) {
      time = (proc.time() - st)['elapsed']
      #diff = results$time[1] - results$time[2]
      cat(sprintf(" Time = %.2f sec\n", time))
    }

    if(!all_agree(results)) {
      p = case$ped
      print(p)
      print(compare(p, 1))

      if(store_bad)
        BAD = c(BAD, list(case))
      else
        break
    }
    # Update progressbar
    if(useProgBar) setTxtProgressBar(pb, i)
  }

  # Close progressbar
  if(useProgBar) close(pb)

  if(store_bad) {
    print(summaryCases(BAD))
    return(BAD)
  }
  invisible(p)
}

summaryCases = function(x) {
  do.call(rbind, lapply(x, function(b)
    as.data.frame(b[c("pedname", "ids", "nall", "geno", "Xchrom", "swapsex",
                      "alsChar", "alsPerm", "alleles", "muttype", "models",
                      "mutrate", "relabel", "reorder")])))
}
