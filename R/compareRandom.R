
#' @importFrom stats runif rpois
randomTestCase = function(ped = NULL, pedname = NULL, ids = NULL, swapsex = NULL,
                          Xchrom = NULL, nall = NULL, alsChar = NULL, alsPerm = NULL,
                          muttype = NULL, relab = NULL, reorder = NULL) {

  if(!requireNamespace("forrel", quietly = TRUE))
    stop2("Package `forrel` is not installed")
  if(!requireNamespace("pedmut", quietly = TRUE))
    stop2("Package `pedmut` is not installed")


  # Special pedigrees
  PEDS = list(SING = singleton(1),
              NUC = nuclearPed(1),
              LIN = linearPed(2),
              HS = halfSibPed(),
              ANC = ancestralPed(2),
              FSM = fullSibMating(1),
              RAND = "random")

  # Specially interesting subsets of genotyped individuals.
  SUBSETS = list(SING = list(),
                 NUC = list(1,2,3,2:3),
                 LIN = list(1, 3, 5),
                 HS = list(4:5),
                 ANC = list(1,5,6,7,c(1,3), c(1,7), c(3,7), c(1,4,7), 5:7),
                 FSM = list(6, 5:6, c(1,6)),
                 RAND = list())


  ### Pedigree
  if(is.null(ped)) {

    # Pedigree
    if(is.null(pedname))
      pedname = sample(names(PEDS), size = 1)

    if(pedname == "RAND") {
      founders = rpois(1, 2) + 1
      matings = sample(founders + 2:3, size= 1)
      ped = randomPed(founders + matings, founders)
    }
    else
      ped = PEDS[[pedname]]

    # permute genders?
    if(is.null(swapsex))
      swapsex = sample(c(T,F), 1)
    if(swapsex) {
      swap_ids = labels(ped)[sample(c(T,F), pedsize(ped), replace = T)]
      ped = swapSex(ped, swap_ids, verbose = F)
    }
  }
  else {
    stopifnot(is.ped(ped))
    ped$MARKERS = NULL
    for(arg in c("pedname", "swapsex", "relab", "reorder"))
      if(!is.null(get(arg))) stop2("When `ped` is given, `", arg, "` must be NULL")
    swapsex = relab = reorder = FALSE
    pedname = NA
  }

  # Subset of genotyped individuals
  if(is.null(ids)) {
    randomsub = labels(ped)[sample(c(T,F), pedsize(ped), replace = T)]
    subsetList = list(None = character(0),
                      Leaves = leaves(ped),
                      All = labels(ped),
                      Random = randomsub)
    if(!is.null(pedname))
      subsetList = c(subsetList, Fixed = SUBSETS[[pedname]])

    ids = sample(subsetList, size = 1)[[1]]
  }
  else {
    stopifnot(all(ids %in% labels(ped)))
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
  x = reorderPed(x, labels(ped))

  m = x$MARKERS[[1]]
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
compareRandom = function(n = 1, programs = NULL, verbose = F, store_bad = FALSE, ...) {
  if(store_bad)
    BAD = list()

  # Progress bar
  useProgBar = n > 9 && !verbose
  if(useProgBar) pb = txtProgressBar(min = 0, max = n, style = 3)

  for(i in 1:n) {
    case = randomTestCase(...)

    if(is.null(programs))
      use_programs = c("pedprobr", ifelse(case$Xchrom, "merlin", "Familias"))
    else
      use_programs = programs

    if(verbose) {
      mutmod = paste0(if(case$models == "") "-" else case$models, ".")
      chr = if(case$Xchrom) "Xchr" else "Auto"
      cat(sprintf("%d. %4s, %s, %d als, mutmod = %-10s ", i, case$pedname, chr, case$nall, mutmod))
      st = proc.time()
    }

    # Run likelihood calculations
    setTimeLimit(elapsed = 10, transient = TRUE)
    results = tryCatch({
      compare(case$ped, 1, programs = use_programs, verbose = FALSE)
      }, error = function(e) NULL)
    setTimeLimit(elapsed = Inf, transient = TRUE)

    if(verbose) {
      time = (proc.time() - st)['elapsed']
      cat(sprintf("Time = %.2f sec (%s)\n", time, paste0(results$program, collapse="/")))
    }

    if(identical(verbose, 2))
      print(results)

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
  invisible(case)
}

summaryCases = function(x) {
  do.call(rbind, lapply(x, function(b)
    as.data.frame(b[c("pedname", "ids", "nall", "geno", "Xchrom", "swapsex",
                      "alsChar", "alsPerm", "alleles", "muttype", "models",
                      "mutrate", "relabel", "reorder")])))
}


