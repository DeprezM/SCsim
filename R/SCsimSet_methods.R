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

newSCsimSet <- function(nGenes, nCells, nPop, pPop, seed = 75,
                        distribution = "gamma", fstParam = 2, sndParam = 0.75,
                        nbBatch, cellsPerBatch, batchEffect,
                        pDEG, pDE, pDM, pDP, pDC, pUp, pDown,
                        cellMixedDM, mixDM,
                        cellMixedDP, mixDP, popMixDP,
                        trajectory, doublet = 2,
                        distrUpFc = "medium", distrDownFc = "medium",
                        dropoutPct) {

  # Check =====================================================================
  ## Checking the validity of numerical variable
  for (dataname in c("nGenes", "nCells", "nPop", "pPop", "seed")) {
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
  for (dataname in c("nGenes", "nCells", "nPop", "seed")) {
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
                       round(1/nPop,2), " %"))
        pPop = rep(round(1/nPop,2), nPop - 1)
        pPop = c(pPop, 1-sum(pPop))
      }
    } else if (sum(pPop) != 1) {
      stop("sum of popProp must be equal to 1")
    }
  }

  # Simulation ================================================================
  # Sample Info
  sampleInfo <- list(nGenes = nGenes,
                     nCells = nCells,
                     nPop = nPop,
                     pPop = pPop)

  # Basal gene expression
  message("-> Compute basal gene expression")
  baseMeans <- computeBaseMean(distribution, fstParam, sndParam,
                               size = nGenes, seed)

  # Batch Effect
  if (!is.null(nbBatch) & nbBatch > 1){
    message("-> Compute batch effect")
    batch <- computeBatch(nbBatch, cellsPerBatch, batchEffect, nCells, seed)
  }

  # Library Size
  message("-> Compute cell biais (library size variation)")
  libSize <- computeLibrarySize(nPop, pPop, nCells, seed)

  # Differentially expressed genes
  message("-> Compute DEG")
  DEG_table <- describeDEG(nGenes, nCells, nPop, pPop,
                           pDEG, pDE, pDP, pDM, pDC,
                           pUp, pDown)

  cell_table <- describeCells(nCells, nPop, pPop, seed,
                            cellMixedDP, cellMixedDM,
                            popMixDP, mixDP, mixDM,
                            trajectory, doublet, DEG_table)

  FC_table <- computeFC(DEG_table, seed, trajectory,
                        popMixDP, distrUpFc, distrDownFc)


  # Dataset Construction ======================================================
  # Effective Means
  message("-> Compute genes effective mean expression")
  effectiveMeans <- replicate(nCells, baseMeans)

  # Apply Batch effect
  if (!is.null(nbBatch) & nbBatch > 1){
    for (b in 1:nbBatch){
      effectiveMeans[ , unlist(batch[[1]][b])] <-
        effectiveMeans[ , unlist(batch[[1]][b])] *
        as.vector(unlist(batch[[2]][b]))
    }
  }

  # Apply Library Size
  effectiveMeans <- effectiveMeans * as.vector(unlist(libSize))
  colnames(effectiveMeans) <- colnames(cell_table)

  # Apply DEG
  for (g in 1:nrow(effectiveMeans)) {
    for (p in 1:nPop) {
      tmpCellIndex <- colnames(cell_table[g + 1, cell_table[g + 1, ] == p])
      effectiveMeans[g, tmpCellIndex] <-
        effectiveMeans[g, tmpCellIndex] * (2^FC_table[g, p])
    }
  }
  rownames(effectiveMeans) <- rownames(FC_table)

  # Apply doublet
  doubletIndex <- colnames(cell_table[1, !cell_table[1, ] %in% 1:nPop] )
  for (i in doubletIndex) {
    effectiveMeans[ ,i] <- rowMeans(effectiveMeans[,
                              unlist(strsplit(cell_table[1, i], "_"))])
  }

  # Basal counts
  message("-> Compute genes basal counts ")
  set.seed(seed)
  baseCounts <- matrix(rnbinom(nGenes*nCells, size = 1,
                               mu = effectiveMeans), ncol = nCells)
  colnames(baseCounts) <- colnames(effectiveMeans)
  rownames(baseCounts) <- rownames(effectiveMeans)

  # Dropout ===================================================================
  message("-> Add dropout events")
  dropout <- matrix(1, ncol = ncol(baseCounts), nrow = nrow(baseCounts))
  colnames(dropout) <- colnames(baseCounts)
  rownames(dropout) <- rownames(baseCounts)
  effectiveCounts <- baseCounts

  isExpr_mat = baseCounts == 0
  pctDropout_genes <- colSums(isExpr_mat) / nrow(isExpr_mat)
  pctDropout_cells <- rowSums(isExpr_mat) / ncol(isExpr_mat)

  set.seed(seed)
  if(median(pctDropout_genes) < (dropoutPct/100)) {
    # rank cells by total counts / library size
    sumCounts <- colSums(baseCounts)
    sumCounts <- sort(sumCounts, decreasing = T)

    # Assign dropout percentage
    cellsDropout <- sort(rnorm(length(sumCounts),
                               mean = dropoutPct/100, sd = 0.05))
    names(cellsDropout) <- names(sumCounts)

    # Estimate the number of genes to dropout
    nbGenesDropout <- ceiling(cellsDropout * nGenes)

    # measure genes mean expression
    genesExpression <- rowMeans(baseCounts)
    names(genesExpression) <- rownames(baseCounts)
    genesExpression <- sort(genesExpression)

    # Assign dropout probability
    # dropoutProba <- (1/(pi * (1 + x^2))) * (1/max(1/(pi * (1 + x^2))))
    # dropoutProba <- 1 * exp(-1 * log2(genesExpression+1))
    lambda <- 0.1
    dropoutProba <- abs(((lambda * exp(-lambda * genesExpression)) /
                      (max(lambda * exp(-lambda * genesExpression)) + 0.025)))

    for (i in 1:length(dropoutProba)) {
      rdm <- runif(1, min = 0, max = 1)
      if (rdm < 0.2) {
        dropoutProba[[i]] <- dropoutProba[[i]] + runif(1, min = 0, max = 0.3)
      } else if (rdm > 0.8) {
        dropoutProba[[i]] <- dropoutProba[[i]] - runif(1, min = 0, max = 0.2)
      }

      if (dropoutProba[[i]] >= 1) {
        dropoutProba[[i]] <- 1 - (dropoutProba[[i]] - 1)
      } else if (dropoutProba[[i]] < 0) {
        dropoutProba[[i]] <- abs(dropoutProba[[i]])
      }
    }
    # plot(log2(genesExpression+1), dropoutProba)

    names(dropoutProba) <- names(genesExpression)
    pb <- txtProgressBar(min = 0, max = nCells, style = 3)

    for (i in colnames(baseCounts)) {
      nbZero <- sum(baseCounts[,i] == 0)
      rdmGenes <- sample(names(dropoutProba), length(dropoutProba))
      g <- 1
      while(nbZero < nbGenesDropout[[i]]) {
        rdm <- runif(1, min = 0, max = 1)
        if (g > length(rdmGenes)) {
          g <- 1
          rdmGenes <- rdmGenes[dropout[rdmGenes, i] !=  0]
        }
        rdmGene <- rdmGenes[g]
        if (rdm < dropoutProba[[rdmGene]]) {
          dropout[rdmGene, i] <- 0
          effectiveCounts[rdmGene, i] <- 0
          nbZero <- sum(effectiveCounts[, i] == 0)
        }
        g <- g + 1
      }
      setTxtProgressBar(pb, grep(i, colnames(baseCounts))[1])
    }
    close(pb)
  }
  isExpr_mat = effectiveCounts == 0
  pctDropout_genes <- colSums(isExpr_mat) / nrow(isExpr_mat)
  pctDropout_cells <- rowSums(isExpr_mat) / ncol(isExpr_mat)
  dropoutPct <- median(pctDropout_genes)

  # Object construction =======================================================
  message("\n     -- Single Cell dataset simulation completed --")

  # Library size dataframe
  cellPopId <- c()
  for (i in 1:nPop){
    cellPopId <- c(cellPopId, rep(i, length(libSize[[i]])))
  }
  lib_df <- rbind(population = as.integer(cellPopId),
                  librarySize = unlist(libSize))

  # Batch dataFrame
  if (!is.null(nbBatch) & nbBatch > 1){
    batchPop <- c()
    for (i in 1:nbBatch){
      batchPop <- c(batchPop, rep(i, length(unlist(batch[[1]][i]))))
    }
    batch_table <- as.data.frame(cbind(batchCell = unlist(batch[[1]]),
                                       batchId = batchPop,
                                       batchValue = unlist(batch[[2]])))

    batch_table <- batch_table[with(batch_table, order(batchCell)), ]
    batch_table <- cbind(population = as.integer(cellPopId),
                         batch_table)
    cell_table <- rbind(batch = batch_table$batchId, cell_table)
  } else {
    batch_table <- data.frame()
  }

  # DEG dataframes
  DEG_df <- cbind(DEG_table[,2:4], FC_table)
  rownames(DEG_df) <- DEG_table$geneLabel


  ##### Generate new SCsimSet object #####
  scsimset <- new("SCsimSet",
                  sampleInfo = sampleInfo,
                  batch_table = batch_table,
                  lib_df = lib_df,
                  DEG_df = DEG_df,
                  cell_table = cell_table,
                  dropout = dropout,
                  baseMeans = baseMeans,
                  effectiveMeans = effectiveMeans,
                  baseCounts = baseCounts,
                  effectiveCounts = effectiveCounts)

  ## Check object validity
  #validObject(distrset)
  scsimset

}

