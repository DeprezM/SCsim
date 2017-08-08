# Distribution methods

computeBaseMean <- function(distribution,
                            fstParam, sndParam, size, seed)
{
  ## Checking valid argument type
  for (dataname in c("fstParam", "sndParam")) {
    eData <- get(dataname)
    if (length(eData) == 0) {
      stop(paste0("argument ", dataname, " must be of length > 0"))
    } else if (!class(eData) %in% c("numeric", "integer")) {
      stop(paste0("argument ", dataname, " must be numeric"))
    }
  }

  if (!missing(seed)) {
    set.seed(seed)
  }

  if (class(distribution) == "character" &
      any(grepl(distribution, c("gamma", "negative binomial")))){

    # Check parameters value depending on the the wanted distribution
    if (distribution == "gamma") {
      if (fstParam <= 0 | sndParam <= 0) {
        stop("fstParam and sndParam must be > 0")
      } else {
        simulated <- rgamma(size, shape=fstParam, rate=sndParam)
      }

    } else if (distribution == "negative binomial") {
      if (fstParam <= 0 | sndParam <= 0) {
        stop("fstParam and sndParam must be > 0")
      } else {
        simulated <- rnbinom(size, size=fstParam, mu=sndParam)
      }

    }
  } else {
    stop("wrong distribution argument, must be one of ('binomial negative',
         'gamma')")
  }

  return(simulated)
}
