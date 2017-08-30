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
                       round(100/nPop,2), " %"))
        pPop = rep(round(100/nPop,2), nPop - 1)
        pPop = c(pPop, 100-sum(pPop))
      }
    } else if (sum(pPop) != 100) {
      stop("sum of popProp must be equal to 1")
    }
  }

  # Simulation ================================================================
  # Sample Info
  sampleInfo <- list(nGenes = nGenes,
                     nCells = nCells,
                     nPop = nPop,
                     pPop = pPop,
                     baseMean = c(distribution, fstParam, sndParam))

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
    lambda <- 0.4
    dropoutProba <- abs(((lambda * exp(-lambda * genesExpression)) /
                      (max(lambda * exp(-lambda * genesExpression)) + 0.025)))

    for (i in 1:length(dropoutProba)) {
      rdm <- runif(1, min = 0, max = 1)
      if (rdm < 0.2) {
        dropoutProba[[i]] <- dropoutProba[[i]] + runif(1, min = 0, max = 0.2)
      } else if (rdm > 0.8) {
        dropoutProba[[i]] <- dropoutProba[[i]] - runif(1, min = 0, max = 0.1)
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
          rdmGenes <- sample(rdmGenes, length(rdmGenes))
        }
        rdmGene <- rdmGenes[g]
        if (rdm < dropoutProba[[rdmGene]]) {
          dropout[rdmGene, i] <- 0
          effectiveCounts[rdmGene, i] <- 0
          nbZero <- sum(effectiveCounts[, i] == 0)
          rdmGenes <- rdmGenes[rdmGenes != rdmGene]
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
  lib_df <- data.frame(population = as.integer(cellPopId),
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

  scsimset
}

# Graphic functions ===========================================================

setMethod("plotDatasetInfo", "SCsimSet",
   function(object) {
     tmpTibble <- tibble::tibble(pPop = object@sampleInfo$pPop,
                                 nPop = as.factor(1:object@sampleInfo$nPop))
     gene <- paste0("  ", object@sampleInfo$nGenes, " Genes  ")
     cell <- paste0("  ", object@sampleInfo$nCells, " Cells  ")

     ggplot2::ggplot(tmpTibble, ggplot2::aes(x="", y=pPop)) +
       ggplot2::geom_bar(ggplot2::aes(fill = nPop), stat = "identity") +
       ggplot2::scale_y_continuous(labels = function(x){ paste0(x, "%") }) +
       ggplot2::coord_flip() +
       ggplot2::ggtitle(bquote(" " %<-% .(cell) %->% " ")) +
       ggplot2::xlab(bquote(" " %<-% .(gene) %->% " ")) +
       ggplot2::labs(y = "") +
       ggplot2::scale_fill_discrete(name="Cell \nPopulation") +
       # ggplot2::theme_bw() +
       ggplot2::theme(
         plot.title = ggplot2::element_text(hjust = 0.5, size = 13),
         axis.text.y = ggplot2::element_blank(),
         axis.title.y = ggplot2::element_text(size = 12),
         axis.ticks.y = ggplot2::element_blank(),
         axis.title.x = ggplot2::element_text( size = 13),
         axis.text.x = ggplot2::element_text(size = 12),
         legend.background = ggplot2::element_rect(fill = "transparent",
                                                   colour = "grey20"),
         legend.position = "bottom",
         plot.margin = grid::unit(c(1.5,1.5,1.5,1.5), "lines"),
         plot.background = ggplot2::element_rect(fill = "transparent"),
         panel.background = ggplot2::element_rect(fill = "transparent"),
         panel.border = ggplot2::element_blank(),
         panel.grid.major = ggplot2::element_blank(),
         panel.grid.minor = ggplot2::element_blank())

   }
)


setMethod("plotBaseMean", "SCsimSet",
   function(object) {
     tmpTibble <- tibble::tibble(simulated = object@baseMeans)
     bin <- (max(tmpTibble$simulated) - min(tmpTibble$simulated))/100
     ggplot2::ggplot(tmpTibble, ggplot2::aes(x=simulated)) +
       ggplot2::geom_histogram(ggplot2::aes(y = ..density..),
                               colour = "deepskyblue3", fill = "deepskyblue1",
                               alpha=1/3, binwidth = bin) +
       ggplot2::geom_line(ggplot2::aes(x=simulated), stat = 'density',
                          colour = "brown4", size = 0.7) +
       ggplot2::annotate("text", y= Inf, x = Inf, hjust = 1, vjust = 1.5,
                         label = paste0("Distribution = ",
                                        object@sampleInfo$baseMean[1]),
                         colour="grey20", size = 4) +
       ggplot2::xlab("Gene expression") +
       ggplot2::ylab("Density") +
       ggplot2::theme(axis.text = ggplot2::element_text(size=11),
         axis.title = ggplot2::element_text(colour="grey20", size = 12),
         axis.line = ggplot2::element_line(colour = "grey30"),
         legend.background = ggplot2::element_rect(fill = "transparent",
                                                   colour = "grey20"),
         plot.background = ggplot2::element_rect(fill = "transparent"),
         panel.background = ggplot2::element_rect(fill = "transparent"),
         panel.border = ggplot2::element_blank(),
         panel.grid.major = ggplot2::element_blank(),
         panel.grid.minor = ggplot2::element_blank())
  })

setMethod("plotBatchInDataset", "SCsimSet",
   function(object) {
     batch_prop <- round(((as.vector(table(object@batch_table[,
           c("population","batchId")])))/object@sampleInfo$nCells)*100)

     tmpTibble <- tibble::tibble(
       pPop = rep(rep(object@sampleInfo$pPop,
                      length(unique(object@batch_table$batchId))),batch_prop),
       nPop = rep(rep(as.factor(1:object@sampleInfo$nPop),
                      length(unique(object@batch_table$batchId))),batch_prop),
       batchID = rep(as.factor(sort(rep(unique(object@batch_table$batchId),
                                        object@sampleInfo$nPop))),batch_prop))


     ggplot2::ggplot(tmpTibble, ggplot2::aes(x=nPop, color = nPop,
                                                    fill = batchID)) +
       ggplot2::geom_bar(stat = "count", size=0.8) +
       ggplot2::scale_fill_brewer(palette="Greys", name="Batch")+
       ggplot2::scale_colour_discrete(name = "Cell \npopulation") +
       ggplot2::scale_y_continuous(labels = function(x){
         paste0(x, "%") }) +
       ggplot2::xlab("Cell Population") +
       ggplot2::labs(y = "Cell population proportion") +
       ggplot2::theme(
         plot.title = ggplot2::element_text(hjust = 0.5, size = 13),
         axis.title.y = ggplot2::element_text(size = 12),
         axis.ticks.y = ggplot2::element_blank(),
         axis.title.x = ggplot2::element_text(colour = "grey30", size = 12),
         axis.text.x = ggplot2::element_text(size = 12),
         axis.text.y = ggplot2::element_text(size = 12),
         legend.background = ggplot2::element_rect(fill = "transparent",
                                                   colour = "grey20"),
         plot.margin = grid::unit(c(1.5,1.5,1.5,1.5), "lines"),
         plot.background = ggplot2::element_rect(fill = "transparent"),
         panel.background = ggplot2::element_rect(fill = "transparent"),
         panel.border = ggplot2::element_blank(),
         panel.grid.major = ggplot2::element_blank(),
         panel.grid.minor = ggplot2::element_blank())

   })


setMethod("plotBatch", "SCsimSet",
   function(object) {
     tmpTibble <- tibble::as.tibble(object@batch_table[,
                                      c("batchId", "batchValue")])
     tmpTibble$batchId <- as.factor(tmpTibble$batchId)

     ggplot2::ggplot(tmpTibble, ggplot2::aes(x = batchValue, fill = batchId)) +
       ggplot2::geom_density(alpha = 0.5) +
       ggplot2::scale_fill_brewer(palette="Greys", name="Batch") +
       ggplot2::xlab("Batch Effect") +
       ggplot2::ylab("Density") +
       ggplot2::theme(
         axis.line = ggplot2::element_line(colour = "grey30"),
         axis.ticks.y = ggplot2::element_blank(),
         axis.title = ggplot2::element_text(colour = "grey30", size = 13),
         axis.text = ggplot2::element_text(size = 12),
         legend.background = ggplot2::element_rect(fill = "transparent",
                                                   colour = "grey20"),
         plot.margin = grid::unit(c(1.5,1.5,1.5,1.5), "lines"),
         plot.background = ggplot2::element_rect(fill = "transparent"),
         panel.background = ggplot2::element_rect(fill = "transparent"),
         panel.border = ggplot2::element_blank(),
         panel.grid.major = ggplot2::element_blank(),
         panel.grid.minor = ggplot2::element_blank())
   })


setMethod("plotCellBiais", "SCsimSet",
   function(object) {
     tmpTibble <- tibble::as.tibble(object@lib_df)
     tmpTibble$population <- as.factor(tmpTibble$population)

     ggplot2::ggplot(tmpTibble, ggplot2::aes(x = librarySize,
                                             fill = population)) +
       ggplot2::geom_density(alpha = 0.5) +
       ggplot2::scale_fill_discrete(name="Population") +
       ggplot2::xlab("Library Effect") +
       ggplot2::ylab("Density") +
       ggplot2::theme(
         axis.line = ggplot2::element_line(colour = "grey30"),
         axis.ticks.y = ggplot2::element_blank(),
         axis.title = ggplot2::element_text(colour = "grey30", size = 13),
         axis.text = ggplot2::element_text(size = 12),
         legend.background = ggplot2::element_rect(fill = "transparent",
                                                   colour = "grey20"),
         plot.margin = grid::unit(c(1.5,1.5,1.5,1.5), "lines"),
         plot.background = ggplot2::element_rect(fill = "transparent"),
         panel.background = ggplot2::element_rect(fill = "transparent"),
         panel.border = ggplot2::element_blank(),
         panel.grid.major = ggplot2::element_blank(),
         panel.grid.minor = ggplot2::element_blank())
   })


setMethod("plotLibrarySize", "SCsimSet",
   function(object) {
     tmpTibble <- tibble::tibble(librarySize = colSums(object@effectiveCounts),
                    population = as.factor(object@lib_df$population))

     ggplot2::ggplot(tmpTibble, ggplot2::aes(x = librarySize,
                                             fill = population)) +
       ggplot2::geom_density(alpha = 0.5) +
       ggplot2::scale_fill_discrete(name="Population") +
       ggplot2::xlab("Library Size") +
       ggplot2::ylab("Density") +
       ggplot2::theme(
         axis.line = ggplot2::element_line(colour = "grey30"),
         axis.ticks.y = ggplot2::element_blank(),
         axis.title = ggplot2::element_text(colour = "grey30", size = 13),
         axis.text = ggplot2::element_text(size = 12),
         legend.background = ggplot2::element_rect(fill = "transparent",
                                                   colour = "grey20"),
         plot.margin = grid::unit(c(1.5,1.5,1.5,1.5), "lines"),
         plot.background = ggplot2::element_rect(fill = "transparent"),
         panel.background = ggplot2::element_rect(fill = "transparent"),
         panel.border = ggplot2::element_blank(),
         panel.grid.major = ggplot2::element_blank(),
         panel.grid.minor = ggplot2::element_blank())
  })

setMethod("plotDEG", "SCsimSet",
   function(object) {
     tmp_pct <- as.data.frame((
       table(object@DEG_df$population)/nrow(object@DEG_df))*100)
     tmpTibble <- tibble::tibble(Percentage = tmp_pct[, 2],
                             Population = paste0("Pop ", (tmp_pct[, 1])),
                             Type = c("noDEG", rep("DEG", nrow(tmp_pct) - 1)))

     tmp_color <- c("#d98c8c", "lightblue", "lightblue",
                    get_color_hexa(nrow(tmp_pct)-1))
     names(tmp_color) <- c("DEG", "noDEG", NA,
                           paste0("Pop ", 1:(nrow(tmp_pct)-1)))

     tmpTibble$Population[tmpTibble$Population == "Pop 0"] <- NA

     ggplot2::ggplot(tmpTibble, ggplot2::aes(x = "", y = Percentage,
                                             fill = Type)) +
       ggplot2::geom_bar(stat = "identity", color = "grey20") +
       ggplot2::scale_fill_manual(name = "Gene type",values = tmp_color) +
       ggplot2::geom_label(ggplot2::aes(label=Population,
                                        fill=factor(Population,
                    levels = sort(tmpTibble$Population, decreasing = T))),
                           na.rm = T, show.legend = F,
                           position = ggplot2::position_stack(vjust = 0.5)) +
       ggplot2::theme(
         axis.text.x = ggplot2::element_blank(),
         axis.title.x = ggplot2::element_blank(),
         axis.ticks.x = ggplot2::element_blank(),
         axis.title.y = ggplot2::element_text(colour = "grey30", size = 12),
         axis.text.y = ggplot2::element_text(size = 11),
         legend.background = ggplot2::element_rect(fill = "transparent",
                                                   colour = "grey20"),
         plot.margin = grid::unit(c(1,1,1,1), "lines"),
         plot.background = ggplot2::element_rect(fill = "transparent"),
         panel.background = ggplot2::element_rect(fill = "transparent"),
         panel.border = ggplot2::element_blank(),
         panel.grid.major = ggplot2::element_blank(),
         panel.grid.minor = ggplot2::element_blank())

   })

setMethod("plotDEGTypes", "SCsimSet",
   function(object) {
     # Create data frame for plot
     tmpTibble <- tibble::tibble(Gene = c("DE", "DM", "DP", "DC"))
     col_n <- c("Gene", paste0("Pop ",
                               1:(length(unique(object@DEG_df$population))-1)))
     for (i in 1:(length(unique(object@DEG_df$population))-1)) {
       tmp_pct <- as.data.frame((
         table(object@DEG_df$type[object@DEG_df$population == i])/
           nrow(object@DEG_df[object@DEG_df$population == i, ]))*100)
       rownames(tmp_pct) <- tmp_pct$Var1
       tmpTibble <- cbind(tmpTibble, tmp_pct[unlist(tmpTibble[,1]),2])
     }
     colnames(tmpTibble) <- col_n

     tmpTibble <- tidyr::gather(tmpTibble, Population, Percentage,
                                col_n[2:length(col_n)])
     tmpTibble$Gene <- factor(tmpTibble$Gene,
                              levels = c("DE", "DM", "DP", "DC"))

     # create color scale
     tmp_color <- c("#42668a", "#6199d1", "#669999", "#3d8f66",
                    get_color_hexa(length(unique(tmpTibble$Population))))
     names(tmp_color) <- c("DE", "DM", "DP", "DC",
                    paste0("Pop ", 1:length(unique(tmpTibble$Population))))

     # plot
     ggplot2::ggplot(tmpTibble, ggplot2::aes(x = Population, y = Percentage,
                                             fill=Gene)) +
       ggplot2::geom_bar(stat = "identity", color = "grey30") +
       ggplot2::scale_fill_manual(values=tmp_color, name = "Gene type") +
       ggplot2::geom_label(ggplot2::aes(label=Population, fill = Population),
                           na.rm = T, show.legend = F,
                           position = ggplot2::position_fill(vjust = 0.5)) +
       ggplot2::theme(
         axis.text.x = ggplot2::element_blank(),
         axis.title.x = ggplot2::element_blank(),
         axis.ticks.x = ggplot2::element_blank(),
         axis.title.y = ggplot2::element_text(colour = "grey50", size = 12),
         axis.text.y = ggplot2::element_text(size = 11),
         legend.background = ggplot2::element_rect(fill = "transparent",
                                                   colour = "grey20"),
         plot.margin = grid::unit(c(1,1,1,1), "lines"),
         plot.background = ggplot2::element_rect(fill = "transparent"),
         panel.background = ggplot2::element_rect(fill = "transparent"),
         panel.border = ggplot2::element_blank(),
         panel.grid.major = ggplot2::element_blank(),
         panel.grid.minor = ggplot2::element_blank())

   })


setMethod("plotCellsContent", "SCsimSet",
   function(object) {
     tmpDF <- object@cell_table
     tmpDF <- tmpDF[3:nrow(tmpDF),]
     # Change DC value - Common DEG
     tmpDF[tmpDF == -1 ] <- object@sampleInfo$nPop + 1
     # Change doublet value
     tmpDF[, ! object@cell_table[2,] %in% 1:object@sampleInfo$nPop] <- -1

     # Gene Label & color
     label <- rownames(tmpDF)
     for (i in 1:length(rownames(tmpDF))) {
       label[i] <- paste(unlist(strsplit(rownames(tmpDF)[i], "_"))[2:3],
                         collapse = "_")
     }
     nDEGtype <- length(unique(label))
     colorLabel <- RColorBrewer::brewer.pal(n=nDEGtype, name="Paired")
     names(colorLabel) <- unique(label)

     # Cell label & color
     labelPop <- as.vector(object@cell_table[2,])
     labelPop[!labelPop %in% 1:object@sampleInfo$nPop] <- "Doublet"
     labelPop[labelPop %in% 1:object@sampleInfo$nPop] <- paste0("Pop ",
                    labelPop[labelPop %in% 1:object@sampleInfo$nPop])

     if (nPop <= 2) {
       colorLabelPop <- RColorBrewer::brewer.pal(n=3,
                           name="YlOrRd")[1:(object@sampleInfo$nPop + 1)]
     } else {
       colorLabelPop <- RColorBrewer::brewer.pal(n=(object@sampleInfo$nPop + 1),
                                                 name="YlOrRd")
     }
     col_common <- colorLabelPop[length(colorLabelPop)]
     colorLabelPop <- c("dodgerblue4",
                        colorLabelPop[1:(length(colorLabelPop)-1)])

     names(colorLabelPop) <- c("Doublet", paste0("Pop ",
                                                 1:object@sampleInfo$nPop))

     annot_col <- data.frame(Pop = unlist(labelPop))
     rownames(annot_col) <- colnames(tmpDF)
     annot_row <- data.frame(DEG = label)
     rownames(annot_row) <- rownames(tmpDF)
     annot_color <- list(DEG = colorLabel,
                         Pop = colorLabelPop)
     colorht <- c("dodgerblue4", "lightblue",
                  colorLabelPop[2:length(colorLabelPop)], col_common)

     pheatmap::pheatmap(tmpDF, color = colorht, border_color = NA,
                        show_rownames = F, show_colnames = F,
                        cluster_rows = F, cluster_cols = F,
                        annotation_row = annot_row, annotation_col = annot_col,
                        annotation_colors = annot_color,
                        annotation_names_row = F, annotation_names_col = F)

   })


setMethod("plotDEGDensity", "SCsimSet",
   function(object) {
     degType <- as.character(unique(scSim@DEG_df$type))
     degType <- degType[degType != "NDE"]
     tmpTibble <- tibble::tibble()

     for (type in degType) {
       fc <- scSim@DEG_df[
         (scSim@DEG_df$type == type) &
           (scSim@DEG_df$regulation == "Up"),]
       if (nrow(fc) > 1) {
         topExpr <- sort(rowMeans(scSim@effectiveCounts[rownames(fc),]),
                         decreasing = T)
         fc <- fc[names(topExpr[1:5]),]
         fc$max <- rowMeans(fc[4:ncol(fc)])
         fc <- fc[with(fc, order(max, decreasing = T)),]
         fc <- fc[rownames(fc[1,]),]
       }
       if (nrow(tmpTibble) == 0) {
         tmpTibble <- tibble::tibble(
           Count = log2(scSim@effectiveCounts[rownames(fc),] + 1),
           Population = as.character(scSim@lib_df$population),
           Type = rep(type, ncol(scSim@effectiveCounts)))
       } else {
         t <- tibble::tibble(
           Count = log2(scSim@effectiveCounts[rownames(fc),] + 1),
           Population = as.character(scSim@lib_df$population),
           Type = rep(type, ncol(scSim@effectiveCounts)))
         tmpTibble <- rbind(tmpTibble, t)
       }
     }

     # NDE case
     fc <- scSim@DEG_df[
       (scSim@DEG_df$type == "NDE"),]
     fc <- sort(rowMeans(scSim@effectiveCounts[rownames(fc),]),
                decreasing = T)
     if (nrow(tmpTibble) == 0) {
       tmpTibble <- tibble::tibble(
         Count = log2(scSim@effectiveCounts[names(fc[1]),] + 1),
         Population = as.character(scSim@lib_df$population),
         Type = rep("NDE", ncol(scSim@effectiveCounts)))
     } else {
       t <- tibble::tibble(
         Count = log2(scSim@effectiveCounts[names(fc[1]),] + 1),
         Population = as.character(scSim@lib_df$population),
         Type = rep("NDE", ncol(scSim@effectiveCounts)))
       tmpTibble <- rbind(tmpTibble, t)
     }

     ggplot2::ggplot(tmpTibble,
                     ggplot2::aes(x = Count, fill = Population)) +
       ggplot2::geom_density(alpha = 0.5) +
       ggplot2::ylim(0,1) +
       ggplot2::facet_grid(~Type)

   })


setMethod("plotDropoutQC", "SCsimSet",
   function(object) {

     isExpr_mat = object@baseCounts == 0
     pctDropout_cells <- rowSums(isExpr_mat) / ncol(isExpr_mat) * 100
     tmpTibble <- tibble::tibble(
       Cells = pctDropout_cells,
       State = factor(rep(paste0("Before ", round(mean(pctDropout_cells),1),
                                 "%"), length(pctDropout_cells))))

     isExpr_mat = object@effectiveCounts == 0
     pctDropout_cells <- rowSums(isExpr_mat) / ncol(isExpr_mat) * 100
     tmpTibble <- rbind(tmpTibble, tibble::tibble(
       Cells = pctDropout_cells,
       State = factor(rep(paste0("After ", round(mean(pctDropout_cells),1),
                                 "%"), length(pctDropout_cells)))))

     p1 = ggplot2::ggplot(tmpTibble, ggplot2::aes(x=Cells)) +
       ggplot2::geom_density(colour = "deepskyblue3", fill = "deepskyblue1",
                             alpha=1/3) +
       ggplot2::xlim(0,100) +
       ggplot2::theme(axis.ticks = ggplot2::element_line(colour="grey",
                                                         size=0.5)) +
       ggplot2::ggtitle("Dropout per cells") +
       ggplot2::labs(x="Dropout percentage") +
       ggplot2::facet_grid(~State)

     # Genes dropout
     isExpr_mat = object@baseCounts == 0
     pctDropout_genes <- colSums(isExpr_mat) / nrow(isExpr_mat) * 100
     tmpTibble <- tibble::tibble(
       Genes = pctDropout_genes,
       State = factor(rep(paste0("Before ", round(mean(pctDropout_genes),1),
                                 "%"), length(pctDropout_genes))))

     isExpr_mat = object@effectiveCounts == 0
     pctDropout_genes <- colSums(isExpr_mat) / nrow(isExpr_mat) * 100
     tmpTibble <- rbind(tmpTibble, tibble::tibble(
       Genes = pctDropout_genes,
       State = factor(rep(paste0("After ", round(mean(pctDropout_genes),1),
                                 "%"), length(pctDropout_genes)))))

     p2 = ggplot2::ggplot(tmpTibble, ggplot2::aes(x=Genes)) +
       ggplot2::geom_density(colour = "deepskyblue3", fill = "deepskyblue1",
                             alpha=1/3) +
       ggplot2::theme_grey() + ggplot2::xlim(0,100) +
       ggplot2::labs(x="Dropout percentage") +
       ggplot2::theme(axis.ticks=ggplot2::element_line(colour="grey", size=0.5)) +
       ggplot2::ggtitle("Dropout per genes") +
       ggplot2::facet_grid(~State)

     scater::multiplot(p1,p2)

   })

setMethod("plotDropoutProba", "SCsimSet",
   function(object) {
     isExpr_mat = object@baseCounts == 0
     pctDropout_cells <- rowSums(isExpr_mat) / ncol(isExpr_mat) * 100
     tmpTibble <- tibble::tibble(
       Cells = pctDropout_cells,
       State = factor(rep(paste0("Before ", round(mean(pctDropout_cells),1),
                                 "%"), length(pctDropout_cells))),
       Expression = log2(rowMeans(object@baseCounts)+1))

     isExpr_mat = object@effectiveCounts == 0
     pctDropout_cells <- rowSums(isExpr_mat) / ncol(isExpr_mat) * 100
     tmpTibble <- rbind(tmpTibble, tibble::tibble(
       Cells = pctDropout_cells,
       State = factor(rep(paste0("After ", round(mean(pctDropout_cells),1),
                                 "%"), length(pctDropout_cells))),
       Expression = log2(rowMeans(object@effectiveCounts)+1)))

     p1 = ggplot2::ggplot(tmpTibble, ggplot2::aes(x=Expression, y=Cells)) +
       ggplot2::geom_point(colour = "deepskyblue3", fill = "deepskyblue1",
                           alpha=1/3) +
       ggplot2::theme_grey() +
       ggplot2::labs(y="Dropout %", x="Mean gene expression") +
       ggplot2::theme(axis.ticks=ggplot2::element_line(colour="grey", size=0.5),
                      text = ggplot2::element_text(size=15)) +
       ggplot2::ggtitle(paste0("% Dropout per cells \nMean: ",
                               round(median(pctDropout_cells),2))) +
       ggplot2::facet_grid(~State)

     # Genes dropout
     isExpr_mat = object@baseCounts == 0
     pctDropout_genes <- colSums(isExpr_mat) / nrow(isExpr_mat) * 100
     tmpTibble <- tibble::tibble(
       Genes = pctDropout_genes,
       State = factor(rep(paste0("Before ", round(mean(pctDropout_genes),1),
                                 "%"), length(pctDropout_genes))),
       Expression = log2(colSums(object@baseCounts)+1))

     isExpr_mat = object@effectiveCounts == 0
     pctDropout_genes <- colSums(isExpr_mat) / nrow(isExpr_mat) * 100
     tmpTibble <- rbind(tmpTibble, tibble::tibble(
       Genes = pctDropout_genes,
       State = factor(rep(paste0("After ", round(mean(pctDropout_genes),1),
                                 "%"), length(pctDropout_genes))),
       Expression = log2(colSums(object@effectiveCounts)+1)))

     p2 = ggplot2::ggplot(tmpTibble, ggplot2::aes(x=Expression, y=Genes)) +
       ggplot2::geom_point(colour = "deepskyblue3", fill = "deepskyblue1",
                           alpha=1/3) +
       ggplot2::theme_grey() +
       ggplot2::labs(y="Dropout %", x="Total counts") +
       ggplot2::theme(axis.ticks=ggplot2::element_line(colour="grey", size=0.5),
                      text = ggplot2::element_text(size=15)) +
       ggplot2::ggtitle(paste0("% Dropout per genes \nMean: ",
                               round(median(pctDropout_genes),2))) +
       ggplot2::facet_grid(~State)


     scater::multiplot(p1,p2)



   })
