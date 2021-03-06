#-----------------------------------------------------------------------------#

computeBatch <- function(nbBatch, cellsPerBatch, batchEffect, nCells, seed) {

  set.seed(seed)

  for (dataname in c("nbBatch", "cellsPerBatch")) {
    eData <- get(dataname)
    if (length(eData) == 0) {
      stop(paste0("argument ", dataname, " must be of length > 0"))
    } else if (!class(eData) %in% c("numeric", "integer")) {
      stop(paste0("argument ", dataname, " must be numeric"))
    }
  }

  if (nbBatch > round(nCells / 100)){
    nbBatch <- as.integer(nCells / 100)
    warning(paste0("Too many batch in this experiment, reset nbBatch = ",
                   nbBatch))
  }

  if (length(cellsPerBatch) >= nbBatch) {
    cellsPerBatch <- cellsPerBatch[1:nbBatch]
  } else if (length(cellsPerBatch) == 1) {
    message(paste0("All batch will have the same size: ",
                   ceiling(1/nbBatch,2)*100, " %"))
    cellsPerBatch <- rep(ceiling(1/nbBatch,2)*100, nbBatch - 1)
    cellsPerBatch = c(cellsPerBatch, 100-sum(cellsPerBatch))
  } else if (length(cellsPerBatch) < nbBatch) {
    stop("argument cellsPerBatch must be of length equal to the number
         of batch")
  }

  if(any(cellsPerBatch < 5) | any(cellsPerBatch > 95)) {
    stop("a minimum of 5% or a maximum of 95% of total percentage of cells per
          batch is required.")
  }

  if(sum(cellsPerBatch) != 100){
    stop("total sum of the percentage of cells per batch must equal 100%")
  }

  if (is.null(batchEffect)) {
    batchEffect <- c(1)
    while (length(batchEffect) < nbBatch) {
      tmpBatch <- sample(c(runif(nbBatch * 2, min = 0.20, max = 0.80),
                           runif(nbBatch * 2, min = 1.20, max = 3.05)), 1)
      if ((length(batchEffect) > 5) |
          (abs(batchEffect[length(batchEffect)] - tmpBatch) > 0.5)) {
        batchEffect <- c(batchEffect, tmpBatch)
      }
    }
  }

  if(length(batchEffect) > nbBatch){
    batchEffect <- batchEffect[1:nbBatch]
  } else if (length(batchEffect) < nbBatch) {
    stop("argument batchEffect must be of length equal to the number
                of batch")
  }

  cellsID <- seq(1, nCells)
  cellsInBatch <- list()
  shiftInBatch <- list()

  for(i in 1:(nbBatch-1)){
    cellsInBatch[[i]] <- sample(cellsID,
                                ceiling((cellsPerBatch[i]/100) * nCells))
    shiftInBatch[[i]] <- abs(rnorm(ceiling((cellsPerBatch[i]/100) * nCells),
                                     mean=batchEffect[i], sd=0.25))
    cellsID <- cellsID[!cellsID %in% cellsInBatch[[i]]]
  }

  cellsInBatch[[nbBatch]] <- cellsID
  shiftInBatch[[nbBatch]] <- abs(rnorm(length(cellsID),
                                         mean=batchEffect[nbBatch], sd=0.25))

  return(list(cellsInBatch, shiftInBatch))
}

#-----------------------------------------------------------------------------#

computeLibrarySize <- function(nPop, pPop, nCells, seed){

  set.seed(seed)
  libSize <- list()

  for (i in 1:(nPop-1)){
    libSize[[i]] <- abs(rnorm(ceiling((pPop[i]/100) * nCells),
                              mean = 1, sd = 0.25))
  }
  libSize[[nPop]] <- abs(rnorm(nCells - length(unlist(libSize)),
                               mean = 1, sd = 0.25))

  return(libSize)
}
