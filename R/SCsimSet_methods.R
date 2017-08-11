# Methods for the SCsimSet class

# Constructor function =========================================================
#' Create a new SCsimSet object.
#'
#' Based on input parameters, it simulates a scRNA-seq dataset (stage by stage):
#' -> basal gene expression
#' -> cell-to-cell biais (~ batch effect + library size)
#'
#'
#' -> FC (DEG)
#' -> effective means (genes mean expression per cell type)
#' -> basal gene counts
#' -> dropout table
#' -> effective counts (with dropout noise)
#'
#' Sample info
#'@param nGenes Scalar of class \code{"numeric"}, providing the number of genes
#' in our simulated dataset.
#'@param nCells Scalar of class \code{"numeric"}, providing the number of cells
#' in our simulated dataset.
#'@param nPop Scalar of class \code{"numeric"}, providing the number of cell
#' population in our simulated dataset.
#'@param pPop Scalar of class \code{"numeric"}, providing the proportion of
#' each cell population in our dataset (percentage of cells per population).
#'@param seed Scalar of class \code{"numeric"}, providing the random number
#' generator (RNG) state for random number generation.
#'
#' Basal expression
#'@param baseDistr **optional** Character string (one of "gamma",
#' "negative binomial") indicating which distribution to use to generate gene
#' expression base means.
#'@param baseFstParam **optional** Scalar of class \code{"numeric"}, providing
#' the 1st distribution parameter. ("gamma" : shape, "negative binomial" : size)
#'@param baseSndParam **optional** Scalar of class \code{"numeric"}, providing
#' the 2nd distribution parameter. ("gamma" : rate, "negative binomial" : mean)
#'
#'@param nbBatch **optional** Scalar of class \code{"numeric"}, providing the
#' number of batch in the dataset (max. a batch every 100 cells).
#'@param cellsPerBatch **optional** Vector of class \code{"numeric"}, providing
#' the number of cells per batch.
#'@param batchEffect **optional** Vector of class \code{"numeric"}, providing
#' the mean shift in library size for each batch.
#'
#'
#' Cell-to-cell biais
#'@param cellDistr **optional** Character string (one of "gamma",
#' "negative binomial", "normal", or "uniform") indicating which distribution to
#' use.
#'@param cellFstParam Scalar of class \code{"numeric"}, providing the 1st
#' distribution parameter. ("gamma" : shape, "negative binomial" : size,
#' "normal" : mean, "uniform" : lower limit)
#'@param cellSndParam Scalar of class \code{"numeric"}, providing the 2nd
#' distribution parameter. ("gamma" : rate, "negative binomial" : mean,
#' "normal" : standard deviation, "uniform" : higher limit).
#'
#' @return a new SCsimSet object containing:
#' \describe{
#'    \item{\code{sampleInfo}:}{Named list which contains sample information:
#'    total number of genes (nGenes), total number of cells (nCells), total
#'    number of cell population (nCells) and relative size / proportion of cell
#'    populations (pPop).}
#'    \item{\code{baseMeans}:}{Object of class \code{"numeric"}, which contains
#'    a vector with genes mean basal expression.}
#'    \item{\code{cellBiais}:}{Object of class \code{"DistrSet"}, which contains
#'    distribution parameters of cell-to-cell biais (capture efficiency, batch,
#'    library size...)}
#'    \item{\code{FcDeg}:}{Object of class \code{"FCSet"}}, which contains DEG
#'    information and FC estimates.
#'    \item{\code{countDispersion}:}{Scalar of class \code{"numeric"}, providing
#'    NB distribution dispersion parameter.}
#'    \item{\code{dropoutPct}:}{Scalar of class \code{"numeric"}, providing the
#'    sample mean dropout percentage.}
#'    \item{\code{dropout}:}{Data Frame which indicates dropout position in
#'    the count matrix.}
#'    \item{\code{effectiveMeans}:}{Data Frame of the dataset effective means
#'    (base means + cell biais + FC).}
#'    \item{\code{baseCounts}:}{Data Frame of the dataset base counts (effective
#'    means + ND estimates).}
#'    \item{\code{effectiveCounts}:}{Data Frame of the dataset effective counts
#'    (baseCounts + dropout).}
#'}
#'
#' @details
#' Parameters set up
#'
#' Basal expression parameters example range of value
#' Gamma distribution: shape = [1.7 : 3], rate = [0.5 : 1.5]
#' Negative Binomial distribution: size = [2 : 3], mu = [2 : 4] # low expression
#'
#' Batch effect : number of batch, maximum 1 batch per 100 cells
#'
#'
#'
#' Based on :
#'      - the number of genes (nGenes)
#'      - the number of cells (nCells)
#'      - the number of cell populations (nPop)
#'      - the proportion of each cell population (pPop)
#'
#' the construction of the SCsimSet object will generate :
#'      -> a descriptive table of the dataset: sampleInfo (list);
#'
#'
#' Based on :
#'      - basal expression distribution parameters (baseDistr, baseFstParam,
#'  baseSndParam)
#'      - cell-to-cell biais distribution parameters (cellDistr, cellFstParam,
#'  cellSndParam)
#'      - the number of genes (nGenes)
#'      - the number of cells (nCells)
#' the construction of the SCsimSet object will generate:
#'      -> baseMeans (DistrSet object)
#'      -> cellBiais (DistrSet object)
#'
#'
#' Based on :
#'      - DEG parameters
#' the construction of the SCsimSet object will generate a FCSet object and
#' effectiveMeans (dataframe)
#'
#'
#' Based on the final dispersion parameter, the construction of the SCsimSet
#' object will generate baseCounts (dataframe)
#'
#'
#' Based on the dropout percentage parameter, the construction of the SCsimSet
#' object will generate:
#'       -> dropout (dataframe)
#'       -> effective counts
#'
#' @export
#' @examples
#'

#'
#'
#'
#'
#'


newSCsimSet <- function(nGenes, nCells, nPop, pPop, seed,
                        baseDistr, baseFstParam, baseSndParam,
                        cellDistr, cellFstParam, cellSndParam,
                        totalDEG, dc, de, dp, dm, up, mix,
                        cellMixed = "constant",
                        upFcDistr = NULL, upFcFstParam = NULL, upFcSndParam = NULL,
                        downFcDistr = NULL, downFcFstParam = NULL, downFcSndParam = NULL,
                        dispersion, dropoutPct)

{

  # Define default parameters
  # defaults <- list("nGenes" = 10000, "nCells" = 500, "nPop" = 2, "pPop" = 0.5,
  #                  "seed" = 1,
  #                  "baseDistr" = "gamma", "baseFstParam" = 2, "baseSndParam" = 0.7,
  #                  "cellDistr" = "normal", "cellFstParam" = -1, "cellSndParam" = 0.25,
  #                  "totalDEG" = 0.1, "dc" = 0.02, "de" = 0.68,
  #                  "dp" = 0.15, "dm" = 0.15, "up" = 0.5, "mix" = 0.3,
  #                  "cellMixed" = "constant",
  #                  "dispersion" = 0.5, "dropoutPct" = 0.8)


  ## Checking if all arguments were provided
  for (dataname in c("nGenes", "nCells", "nPop", "pPop", "seed",
                     "baseDistr", "baseFstParam", "baseSndParam",
                     "cellDistr", "cellFstParam", "cellSndParam",
                     "totalDEG", "dc", "de", "dp", "dm", "up", "mix",
                     "cellMixed", "dispersion", "dropoutPct")) {
    eData <- get(dataname)
    if (missing(eData) | is.null(eData)) {
      warning(paste0("missing argument : ", dataname))
      warning(paste0("Will use default parameter: ", dataname, " = ",
                     defaults[[dataname]]))
      dataname <- defaults[[dataname]]
    }
  }

  #### Sample information ####
  ## Checking the validity of numerical variable
  for (dataname in c("nGenes", "nCells", "nPop", "pPop")) {
    eData <- get(dataname)
    if (length(eData) == 0) {
      stop(paste0("argument ", dataname, " must be of length > 0"))
    } else if (!class(eData) %in% c("numeric", "integer")) {
      stop(paste0("argument ", dataname, " must be numeric"))
    } else if (any(eData <= 0)){
      stop(paste0("argument ", dataname, "must be > 0"))
    }
  }

  ## Checking the validity of number of element (n) variable
  for (dataname in c("nGenes", "nCells", "nPop")) {
    eData <- get(dataname)
    if (eData %% 1 != 0) {
      stop(paste0("argument ", dataname, " must be a whole number (integer)"))
    }
  }

  ## Check validity of nGenes
  if (nGenes < 1000) {
    stop("Simulated dataset must contain at least 1000 genes")
  }

  ## Check validity of nCells
  if (nCells < 100) {
    stop("Simulated dataset must contain at least 100 cells")
  }

  ## Check validity of nPop
  if (nPop == 1) {
    message("! There is only one cell population in the simulated dataset.")
    pPop = 1
  } else {
    ## Check validity of pPop
    if (length(pPop) != nPop) {
      if (length(pPop) == 1) {
        message(paste0("! All cell population will have the same size: ",
                       round(1/nPop,2)*100, " %"))
        pPop = rep(round(1/nPop,2)*100, nPop - 1)
        pPop = c(pPop, 100-sum(pPop))
      }
    } else if (sum(pPop) != 1) {
      stop("Sum of pProp must be equal to 1")
    }
  }

  ## Create sampleInfo dataframe
  sampleInfo <- list(nGenes = nGenes,
                     nCells = nCells,
                     nPop = nPop,
                     pPop = pPop)


  #### baseMeans ####
  message("-> Compute basal gene expression")
  baseMeans <- computeBaseMean(distribution, fstParam, sndParam, size, seed)

  message("-> Compute batch effect")
  baseMeans <- computeBatch(nbBatch, cellsPerBatch, batchEffect, nCells, seed)

  message("-> Compute cell biais (library size variation)")
  baseMeans <- computeLibrarySize(nPop, pPop, nCells, seed)

  message("-> Compute DEG")
