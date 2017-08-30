# describeDEG =================================================================

#' describeDEG -> Differentially expressed genes indexes / type / regulation in
#' the gene expression table to generate
#'

# DEG_table <- describeDEG(nGenes = 1000, nCells = 100, nPop = 3,
#                          pPop = c(15, 35, 50), pDEG = c(20, 30, 10),
#                          pDE = c(40, 10, 20), pDP = 20, pDM = c(30, 60, 50), pDC = 10,
#                          pUp = c(70, 50, 20), pDown= c(30, 50, 80))


describeDEG <- function(nGenes, nCells, nPop, pPop,
                        pDEG, pDE = 0, pDP = 0, pDM = 0, pDC = 0,
                        pUp = 50, pDown= 50) {

  ## Check pDEG
  if (missing(pDEG) | is.null(pDEG)){
    stop("missing argument : pDEG")
  } else if (any(pDEG == 0)) {
    message("no DEG in this dataset")
  }

  if (any(pDEG != 0) | nPop != 1) {
    # Check parameter type
    for (dataname in c("pDE", "pDP", "pDM", "pDC", "pUp", "pDown")){
      eData <- get(dataname)
      if (class(eData) != "numeric" & class(eData) != "integer"){
        stop(paste0("argument : ", dataname, " must be of class numeric."))
      } else if (any(eData < 0) | any(eData > 100)){
        stop(paste0("argument : ", dataname, " is a percentage,
                    which range is between 0 and 100."))
      }
    }

    # Check percentage validity
    if ((length(pDE) != length(pDM))) {
      stop(paste0("pDE and pDM must be of equal length either 1 or ", nPop))
    }

    if ((length(pUp) != length(pDown))) {
      stop(paste0("pUp and pDown must be of equal length either 1 or ", nPop))
    }

    if (length(pDE) > 1 & length(pDM) > 1){
      for (i in 1:nPop){
        if (!isTRUE(all.equal(pDC + pDE[i] + pDP + pDM[i], 100))){
        stop("sum of pDC, pDE, pDP and  pDM != 100. Those proportions are
             sub-types of DEG which sum must be equal to 100.")
        }
      }
    } else {
      if (!isTRUE(all.equal(pDC + pDE + pDP + pDM, 100))){
      stop("sum of pDC, pDE, pDP and  pDM != 100. Those proportions are
             sub-types of DEG which sum must be equal to 100.")
      } else {
        pDE <- rep(pDE, nPop)
        pDM <- rep(pDM, nPop)
      }
    }

    # Check percentage validity
    if (length(pUp) > 1) {
      for(i in 1:nPop) {
        if (!isTRUE(all.equal(pUp[i] + pDown[i], 100))) {
          stop("sum of pUp and  pDown != 100.")
        }
      }
    } else {
      if (!isTRUE(all.equal(pUp + pDown, 100))) {
        stop("sum of pUp and  pDown != 100.")
      } else {
        pUp <- rep(pUp, nPop)
        pDown <- rep(pDown, nPop)
      }
    }

    pUp <- pUp / 100
    pDown <- pDown / 100

    if (length(pDEG) > 1 & length(pDEG) == nPop){
      nDEGpop <- (pDEG/100) * nGenes
    } else if (length(pDEG) == 1) {
      nDEGpop <- rep((pDEG/100) * nGenes, npop)
    } else {
      stop(paste0("pDEG must be of length 1 or ", nPop, " (nPop)"))
    }

    listDEGtype <- c()
    for (dataname in c("pDE", "pDP", "pDM", "pDC")){
      eData <- get(dataname)
      if (any(eData != 0)){
        listDEGtype <- c(listDEGtype, dataname)
      }
    }

    genesStatus <- data.frame(
      index = seq(1, nGenes),
      type = factor(rep("NDE", nGenes),
                    levels=c("NDE", "DC", "DE", "DP", "DM")),
      regulation = factor(rep("none", nGenes), levels=c("Up", "Down", "none")),
      population = rep(0, nGenes),
      geneLabel = rep("", nGenes),
      stringsAsFactors = F
    )

    popIndex <- 1
    for (pop in 1:nPop) {
      genesStatus$population[popIndex:(popIndex + nDEGpop[pop] - 1)] <-
        rep(pop, nDEGpop[pop])
      geneType <- c()
      geneReg <- c()

      if ("pDE" %in% listDEGtype) {
        if (pop == 1) {
          pDE <- pDE/100
        }
        geneType <- c(geneType, rep("DE", ceiling(nDEGpop[pop] * pDE[pop])))
        geneReg <- c(geneReg, rep("Up",
                        ceiling(ceiling(nDEGpop[pop] * pDE[pop]) * pUp[pop])))
        geneReg <- c(geneReg, rep("Down", sum(geneType == "DE") -
                        ceiling(ceiling(nDEGpop[pop] * pDE[pop]) * pUp[pop])))
      }

      if ("pDM" %in% listDEGtype) {
        if (pop == 1) {
          pDM <- pDM/100
        }
        geneType <- c(geneType, rep("DM", ceiling(nDEGpop[pop] * pDM[pop])))
        geneReg <- c(geneReg, rep("Up",
                        ceiling(ceiling(nDEGpop[pop] * pDM[pop]) * pUp[pop])))
        geneReg <- c(geneReg, rep("Down", sum(geneType == "DM") -
                        ceiling(ceiling(nDEGpop[pop] * pDM[pop]) * pUp[pop])))

      }

      if ("pDP" %in% listDEGtype) {
        if (pop == 1) {
          pDP <- pDP/100
        }
        geneType <- c(geneType, rep("DP", ceiling(nDEGpop[pop] * pDP)))
        geneReg <- c(geneReg, rep("Up",
                         ceiling(ceiling(nDEGpop[pop] * pDP) * pUp[pop])))
        geneReg <- c(geneReg, rep("Down", sum(geneType == "DP") -
                         ceiling(ceiling(nDEGpop[pop] * pDP) * pUp[pop])))
      }
      if ("pDC" %in% listDEGtype) {
        if (pop == 1) {
          pDC <- pDC/100
        }
        geneType <- c(geneType, rep("DC", ceiling(nDEGpop[pop] * pDC)))
        geneReg <- c(geneReg, rep("Up",
                         ceiling(ceiling(nDEGpop[pop] * pDC) * pUp[pop])))
        geneReg <- c(geneReg, rep("Down", sum(geneType == "DC") -
                         ceiling(ceiling(nDEGpop[pop] * pDC) * pUp[pop])))
      }
      geneType <- c(rep("DE", nDEGpop[pop]-length(geneType)), geneType)
      geneReg <- c(rep("Up", nDEGpop[pop]-length(geneType)), geneReg)

      genesStatus$type[popIndex:(popIndex + nDEGpop[pop] - 1)] <- geneType
      genesStatus$regulation[popIndex:(popIndex + nDEGpop[pop] - 1)] <- geneReg

      popIndex <- popIndex + nDEGpop[pop]
    }

    for (g in 1:nGenes) {
      genesStatus$geneLabel[g] = paste(genesStatus$index[g],
                                       genesStatus$type[g],
                                       genesStatus$regulation[g],
                                       genesStatus$population[g], sep="_")
    }

  } else {
    genesStatus <- data.frame(
      index = seq(1, nGenes),
      type = factor(rep("NDE", nGenes),
                    levels=c("NDE", "DC", "DE", "DP", "DM")),
      regulation = factor(rep("none", nGenes), levels=c("Up", "Down", "none")),
      population = rep(0, nGenes),
      geneLabel = rep("", nGenes),
      stringsAsFactors = F
    )

    for (g in 1:nGenes) {
      genesStatus$geneLabel[g] = paste(genesStatus$index[g],
                                       genesStatus$type[g],
                                       genesStatus$regulation[g],
                                       genesStatus$population[g], sep="_")
    }

  }

  return(genesStatus)
}


# describeCells ===============================================================

#' describeCells -> Associate cells with their DEG.
# cell_table <- describeCells(nCells = 100, nPop = 3, pPop = c(15,35,50),
#                           cellMixedDP = "pseudo", cellMixedDM = "pseudo",
#                           seed = NULL,
#                           popMixDP = NULL, mixDP = 25, mixDM = 25,
#                           trajectory = list(c(1,2,3)), doublet = 0, DEG_table)


describeCells <- function(nCells, nPop, pPop, seed,
                          cellMixedDP = "pseudo", cellMixedDM = "pseudo",
                          popMixDP = NULL, mixDP = 25, mixDM = 25,
                          trajectory = NULL, doublet = 2, genesStatus) {

  ## Checking arguments presence
  if (is.null(genesStatus)) {
    stop("missing argument : genesStatus")
  }

  ## Checking the validity of numerical variable
  for (dataname in c("mixDP", "mixDM", "doublet")) {
    eData <- get(dataname)
    if (length(eData) == 0) {
      stop(paste0("argument ", dataname, " must be of length > 0"))
    } else if (!class(eData) %in% c("numeric", "integer")) {
      stop(paste0("argument ", dataname, " must be numeric"))
    } else if (any(eData < 0) | (any(eData > 100))){
      stop(paste0("argument ", dataname, "must in range 0 to 100 (percentage)"))
    }
  }

  ## Checking the validity of number of element (n) variable
  for (dataname in c("cellMixedDP", "cellMixedDM")) {
    eData <- get(dataname)
    if (class(eData) != "character") {
      stop(paste0("argument ", dataname, " must be of class character"))
    }
    if (!eData %in% c("random", "constant", "pseudo")) {
      stop(paste0("argument ", dataname, " must be one of pseudo,
                  constant or random"))
    }
  }

  ## Check seed
  if (!is.null(seed) & !class(seed) %in% c("numeric", "integer")){
    stop("seed argument must be integer")
  } else {
    seed <- 75
  }

  ## Check validity of popMixDP
  if (nPop > 2 & !is.null(popMixDP)) {
    nbDP <- length(unique(genesStatus[genesStatus$type == "DP", "population"]))

    if (class(popMixDP) != "list" | length(popMixDP) != nDP) {
      stop(paste0("popMixDP must be a list of length ", nDP))
    }
  }

  if (any(!popMixDP %in% 1:nPop)) {
    stop(paste0("popMixDP must only include population ID in range 1 to ",
                nPop))
  }

  for (i in 1:nPop) {
    if (!is.null(popMixDP))
      if ((length(popMixDP[[i]]) < 2 |
                              length(popMixDP[[i]]) > nPop - 1)) {
      stop(paste0("element ", i, " of popMixDP must be in range from 2 to ",
                  nPop - 1))
    }
  }

  # Check trajectory
  if (nPop > 2 & !is.null(trajectory)) {
    if (class(trajectory)!= "list" |
        (length(trajectory) < 1 | length(trajectory) > 3)) {
      stop("popMixDP must be a list of maximum length 3")
    } else {
      for (i in 1:length(trajectory)) {
        if (length(trajectory[[i]]) < 3 | length(trajectory[[i]]) > nPop) {
          stop(paste0("trajectory must be of length at least 3 to ", nPop))
        }
      }

      for (i in 1:length(trajectory)) {
        if (length(trajectory[[i]]) != length(unique(trajectory[[i]]))) {
          stop("a trajectory cannot have 2 identical transitionnary state")
        }
      }
    }
  }

  if (!is.null(trajectory) & is.null(popMixDP)) {
    tmp_popMixedDP <- list()
    for (i in 1:length(trajectory)){
      for (j in 1:(length(trajectory[[i]])-1)){
        tmp_popMixedDP[[length(tmp_popMixedDP)+1]] <-
          c(trajectory[[i]][j], trajectory[[i]][j+1])
      }
      # Add mixing for last step of trajectory
      tmp_popMixedDP[[length(tmp_popMixedDP)+1]] <-
        c(trajectory[[i]][length(trajectory[[i]])],
          trajectory[[i]][length(trajectory[[i]])-1])
    }
    length(tmp_popMixedDP)
    tmp_popMixedDP <- unique(tmp_popMixedDP)
    starts <- unlist(tmp_popMixedDP)[seq(1,length(tmp_popMixedDP)*2,2)]
    popMixDP <- tmp_popMixedDP[!duplicated(starts)]
    tmp_popMixedDP <- tmp_popMixedDP[duplicated(starts)]

    if (length(tmp_popMixedDP) != 0){
      starts <- unlist(tmp_popMixedDP)[seq(1,length(tmp_popMixedDP)*2,2)]
      for (i in 1:length(popMixDP)){
        if (popMixDP[[i]][1] %in% starts){
          popMixDP[[i]] <- c(popMixDP[[i]],
                             tmp_popMixedDP[[grep(popMixDP[[i]][1], starts)]][2])
        }
      }
    }
  }

  starts <- c()
  for (i in 1:length(popMixDP)){
    starts <- c(starts, popMixDP[[i]][1])
  }
  pPop <- pPop/100

  nGenes <- nrow(genesStatus)

  ## Check validity of nPop
  if (nPop == 1) {
    message("There is only one cell population in simulated dataset,")

    ## Create/Initiate cellsStatus table
    cellsStatus <- as.data.frame(matrix(data = 1,
                                        nrow = nGenes + 1,
                                        ncol = nCells))
    rownames(cellsStatus) = c("cellPop", genesStatus$geneLabel)

  } else {

    ## Create/Initiate cellsStatus table
    cellsStatus <- as.data.frame(matrix(data = 0,
                                        nrow = nGenes + 1,
                                        ncol = nCells))
    rownames(cellsStatus) = c("cellPop", genesStatus$geneLabel)

    # Label each cell with its population number
    tmpIndexCell <- 1
    tmpCntCell <- 0 # Sum of the cell in each population (check)
    for (x in 1:nPop) {
      if (x < nPop) {
        tmpNbCells <- round(nCells*pPop[x])
        tmpCntCell <- tmpCntCell + tmpNbCells
      } else {
        tmpNbCells <- nCells - tmpCntCell
        tmpCntCell <- tmpCntCell + tmpNbCells
      }
      cellsStatus["cellPop", tmpIndexCell:(tmpIndexCell + tmpNbCells - 1)] <-
        rep(x, length(tmpIndexCell:(tmpIndexCell + tmpNbCells - 1)))
      tmpIndexCell <- tmpIndexCell + tmpNbCells
    }

    # Associate each cell with its DEG
    for (i in 1:nPop){
      nDEGtype <- as.character(unique(genesStatus[
        genesStatus$population == i, "type"]))
      set.seed(seed)
      for (type in nDEGtype){
        if (type == "DC") {
          geneName <- genesStatus[(genesStatus$population == i) &
                                    genesStatus$type == type, "geneLabel"]
          cellsStatus[geneName, cellsStatus["cellPop", ] == i] <- -1
        }

        if (type == "DE") {
          geneName <- genesStatus[(genesStatus$population == i) &
                                    genesStatus$type == type, "geneLabel"]
          cellsStatus[geneName, cellsStatus["cellPop", ] == i] <- i
        }

        if (type == "DM") {
          geneName <- genesStatus[(genesStatus$population == i) &
                                    genesStatus$type == type, "geneLabel"]

          if (cellMixedDM == "constant"){
            cell <- sample(colnames(
              cellsStatus[, cellsStatus["cellPop", ] == i]),
              round((mixDM / 100) * sum(cellsStatus["cellPop", ] == i)))
            cellsStatus[geneName, cellPop[!cellPop %in% cell] ] <- i

          } else if (cellMixedDM == "random"){
            for (g in geneName){
              cell <- sample(colnames(
                cellsStatus[, cellsStatus["cellPop", ] == i]),
                round((mixDM / 100) * sum(cellsStatus["cellPop", ] == i)))
              cellsStatus[g, cellPop[!cellPop %in% cell] ] <- i
            }

          } else if (cellMixedDM == "pseudo"){
            fixedCell <- sample(colnames(
              cellsStatus[, cellsStatus["cellPop", ] == i]),
              round((mixDM / 100) * sum(cellsStatus["cellPop", ] == i)))

            cellPop <- colnames(cellsStatus[, cellsStatus["cellPop", ] == i])

            for (g in geneName){
              celltmp <- sample(cellPop[!cellPop %in% fixedCell],
                ceiling(length(fixedCell) / 2))
              cell <- c(celltmp, sample(fixedCell,
                    floor(length(fixedCell) / 2)))
              cellsStatus[g, cellPop[!cellPop %in% cell]] <- i
            }
          }
        }

        if (type == "DP") {
          geneName <- genesStatus[(genesStatus$population == i) &
                                    genesStatus$type == type, "geneLabel"]

          if (i %in% starts){
            mixCell <- popMixDP[starts == i][[1]][2:length(
              popMixDP[starts == i][[1]])]

            if (cellMixedDP == "constant"){
              cell <- sample(colnames(
                cellsStatus[, cellsStatus["cellPop", ] == i]),
                round(((100 - mixDP) / 100) * sum(cellsStatus["cellPop", ] == i)))
              for(c in mixCell){
                cell <- c(cell, sample(colnames(
                  cellsStatus[, cellsStatus["cellPop", ] == c]),
                  round((mixDP / 100) * sum(cellsStatus["cellPop", ] == c))))
              }
              cellsStatus[geneName, cell ] <- i

            } else if (cellMixedDP == "random"){
              for (g in geneName){
                cell <- sample(colnames(
                  cellsStatus[, cellsStatus["cellPop", ] == i]),
                  round(((100 - mixDP) / 100) * sum(cellsStatus["cellPop", ] == i)))
                for(c in mixCell){
                  cell <- c(cell, sample(colnames(
                    cellsStatus[, cellsStatus["cellPop", ] == c]),
                    round((mixDP / 100) * sum(cellsStatus["cellPop", ] == c))))
                }
                cellsStatus[geneName, cell ] <- i
              }

            } else if (cellMixedDP == "pseudo"){
              # cell in given population that will have mix values
              fixedCell <- list(sample(colnames(
                cellsStatus[, cellsStatus["cellPop", ] == i]),
                round(((100 - mixDP)/ 100) * sum(cellsStatus["cellPop", ] == i))))

              cellPop <- colnames(cellsStatus[, cellsStatus["cellPop", ] == i])
              nbPopCell <- round((mixDP / 100) * sum(
                cellsStatus["cellPop", ] == i))

              # cell in related population that will have mix values
              for(c in mixCell){
                fixedCell <- c(fixedCell, list(sample(colnames(
                  cellsStatus[, cellsStatus["cellPop", ] == c]),
                  round((mixDP / 100) * sum(cellsStatus["cellPop", ] == c)))))
              }

              # for each DP genes
              for (g in geneName){
                # fixed cells mixing
                cell <- sample(cellPop[cellPop %in% fixedCell[[1]]],
                               ceiling(nbPopCell/2))

                # random cells mixing
                cell <- c(cell, sample(cellPop[!cellPop %in% fixedCell],
                                       floor(nbPopCell/2)))

                for (c in mixCell){
                  tmpcellPop <- colnames(cellsStatus[, cellsStatus["cellPop",
                                                                   ] == c])
                  # fixed cells mixing
                  cell <- c(cell, sample(tmpcellPop[tmpcellPop %in%
                                     fixedCell[[1+grep(c, mixCell)]]],
                      round((length(fixedCell[[1+grep(c, mixCell)]])/2))))
                  cell <- c(cell, sample(tmpcellPop[!tmpcellPop %in%
                                     fixedCell[[1+grep(c, mixCell)]]],
                      round((length(fixedCell[[1+grep(c, mixCell)]])/2))))
                }
                cellsStatus[g, cellPop[!cellPop %in% cell ]] <- i
                cellsStatus[g, cell[!cell %in% cellPop ]] <- i

              } # for geneName
            } # if cellMix
          } # if i in starts
        } # if geneType
      } # for type in DEG
    } # for cell population
    if (doublet != 0){
      dbCell <- sample(colnames(cellsStatus), floor((doublet / 100)*nCells))
      tmpDb <- sample(colnames(cellsStatus)[!colnames(cellsStatus) %in% dbCell],
                      floor((doublet / 100)*nCells)*3)
      for (c in dbCell){
        pick <- sample(tmpDb, round(runif(1, min = 2, max = 3)))
        cellsStatus["cellPop", c] <- paste(pick, collapse = "_")
        tmpDb <- tmpDb[!tmpDb %in% pick]
      }
    }
  } # if nPop
  return(cellsStatus)
}

# computeFC ===================================================================

# FC_table <- computeFC(genesStatus, seed = 25, trajectory = list(c(1,2,3)),
#                       popMixDP = NULL,
#                       distrUpFc = "medium", distrDownFc = "medium")

computeFC <- function(genesStatus, seed, trajectory = NULL, popMixDP = NULL,
                      distrUpFc = "medium", distrDownFc = "medium") {

  ## Checking arguments presence
  if (is.null(genesStatus)) {
    stop("missing argument : genesStatus")
  }

  ## Get simulated dataset dimension
  nPop <- max(genesStatus$population)
  nGenes <- nrow(genesStatus)
  nDegUp <- sum(as.character(genesStatus$regulation) == "Up")
  nDegDown <- sum(as.character(genesStatus$regulation) == "Down")

  ## Check validity of distr.. arguments
  if (!distrUpFc %in% c("small", "medium", "high")){
    stop("argument distrUpFc must be one of small / medium / high")
  }

  if (!distrDownFc %in% c("small", "medium", "high")){
    stop("argument distrDownFc must be one of small / medium / high")
  }

  set.seed(seed)

  if (distrUpFc == "small"){
    UpFc <- rnbinom(mu = 8, size = 20, n = nDegUp) / 10
  } else if (distrUpFc == "medium") {
    UpFc <- rnbinom(mu = 12, size = 24, n = nDegUp) / 10
  } else if (distrUpFc == "high") {
    UpFc <- rnbinom(mu = 16, size = 28, n = nDegUp) / 10
  }

  if (distrDownFc == "small"){
    DownFc <- -(rnbinom(mu = 8, size = 22, n = nDegDown) / 10)
  } else if (distrDownFc == "medium") {
    DownFc <- -(rnbinom(mu = 12, size = 26, n = nDegDown) / 10)
  } else if (distrDownFc == "high") {
    DownFc <- -(rnbinom(mu = 16, size = 30, n = nDegDown) / 10)
  }


  ## Create/ Initialize trueFc dataframe
  trueFc <- data.frame(matrix(data = 1,
                              nrow = nGenes,
                              ncol = nPop))
  rownames(trueFc) <- genesStatus$geneLabel

  # Complete trueFc dataframe for DE, DP and DM DEG type
  set.seed(seed)
  for (x in 1:nPop) {
    tmpIndexDEGup <- genesStatus[genesStatus$population == x &
                          genesStatus$regulation == "Up", "geneLabel"]
    tmpIndexDEGdown <- genesStatus[genesStatus$population == x &
                            genesStatus$regulation == "Down", "geneLabel"]
    trueFc[tmpIndexDEGup, x] <- sample(UpFc, length(tmpIndexDEGup))
    trueFc[tmpIndexDEGdown, x] <- sample(DownFc, length(tmpIndexDEGdown))
  }

  # Complete trueFc dataframe for DC DEG type
  if (any(grepl("DC", unique(as.character(genesStatus$type))))) {
    if (!is.null(trajectory)){
      for (x in 1:nPop){
        nbTraj <- grep(x, trajectory)
        tmpDC <- genesStatus[genesStatus$population == x &
                               genesStatus$type == "DC", "geneLabel"]
        if (length(nbTraj) > 1){
          trajIndex <- list()
          for (i in 1:(length(nbTraj) - 1)) {
            trajIndex <- c(trajIndex, list(sample(tmpDC,
                               floor(length(tmpDC)/length(nbTraj)))))
            tmpDC <- tmpDC[!tmpDC %in% trajIndex[[i]]]
          }
          trajIndex <- c(trajIndex, list(tmpDC))
        } else {
          trajIndex <- list(tmpDC)
        }

        for (i in 1:length(nbTraj)) {
          tmpTraj <- trajectory[[nbTraj[i]]]
          # Increase expression
          trajIndexInc <- trajIndex[[i]][1:floor(length(trajIndex[[i]])/2)]
          # Decrease expression
          trajIndexDec <- trajIndex[[i]][floor(length(trajIndex[[i]])/2):
                                           length(trajIndex[[i]])]

          # Down half
          for (t in trajIndexInc) {
            tmp_sample <- sample(DownFc,
                                 length(tmpTraj[1:floor(length(tmpTraj)/2)]))
            if (any(duplicated(tmp_sample))) {
              tmp_sample[duplicated(tmp_sample)] <-
                2*tmp_sample[duplicated(tmp_sample)]
            }
            trueFc[t, tmpTraj[1:floor(length(tmpTraj)/2)]] <- sort(tmp_sample)

            tmp_sample <- sample(UpFc, length(tmpTraj[(ceiling(
              length(tmpTraj)/2)+1):length(tmpTraj)]))
            if (any(duplicated(tmp_sample))) {
              tmp_sample[duplicated(tmp_sample)] <-
                2*tmp_sample[duplicated(tmp_sample)]
            }
            trueFc[t, tmpTraj[(ceiling(length(tmpTraj)/2)
                               +1):length(tmpTraj)]] <- sort(tmp_sample)
          }

          for (t in trajIndexDec) {
            tmp_sample <- sample(UpFc,
                                 length(tmpTraj[1:floor(length(tmpTraj)/2)]))
            if (any(duplicated(tmp_sample))) {
              tmp_sample[duplicated(tmp_sample)] <-
                2*tmp_sample[duplicated(tmp_sample)]
            }
            if (length(tmp_sample) != 1) {
              trueFc[t, tmpTraj[1:floor(length(tmpTraj)/2)]] <- sort(tmp_sample,
                                                                decreasing = T)
            } else {
              trueFc[t, tmpTraj[1:floor(length(tmpTraj)/2)]] <- tmp_sample
            }
            tmp_sample <- sample(DownFc, length(tmpTraj[(
              ceiling(length(tmpTraj)/2)+1):length(tmpTraj)]))
            if (any(duplicated(tmp_sample))) {
              tmp_sample[duplicated(tmp_sample)] <-
                2*tmp_sample[duplicated(tmp_sample)]
            }
            if (length(tmp_sample) != 1) {
              trueFc[t, tmpTraj[(ceiling(length(tmpTraj)/2)
                     +1):length(tmpTraj)]] <- sort(tmp_sample, decreasing = T)
            } else {
              trueFc[t, tmpTraj[(ceiling(length(tmpTraj)/2)
                                 +1):length(tmpTraj)]] <- tmp_sample
            }


          }
        }
      }
    } else if (!is.null(popMixDP)) {
      for (x in 1:nPop){
        nbTraj <- grep(x, popMixDP)
        tmpDC <- genesStatus[genesStatus$population == x &
                               genesStatus$type == "DC", "geneLabel"]
        if (length(nbTraj) > 1){
          trajIndex <- list()
          for (i in 1:(length(nbTraj) - 1)) {
            trajIndex <- c(trajIndex, list(sample(tmpDC,
                              floor(length(tmpDC)/length(nbTraj)))))
            tmpDC <- tmpDC[!tmpDC %in% trajIndex[[i]]]
          }
          trajIndex <- c(trajIndex, list(tmpDC))
        } else {
          trajIndex <- list(tmpDC)
        }

        for (i in 1:length(nbTraj)) {
          tmpTraj <- popMixDP[[nbTraj[i]]]
          # Increase expression
          trajIndexInc <- trajIndex[[i]][1:floor(length(trajIndex[[i]])/2)]
          # Decrease expression
          trajIndexDec <- trajIndex[[i]][floor(length(trajIndex[[i]])/2):
                                           length(trajIndex[[i]])]

          # Down half
          for (t in trajIndexInc) {
            tmp_sample <- sample(UpFc, length(tmpTraj))
            if (any(duplicated(tmp_sample))) {
              tmp_sample[duplicated(tmp_sample)] <-
                2*tmp_sample[duplicated(tmp_sample)]
            }
            trueFc[t, tmpTraj] <- sort(tmp_sample)
          }

          for (t in trajIndexDec) {
            tmp_sample <- sample(DownFc, length(tmpTraj))
            if (any(duplicated(tmp_sample))) {
              tmp_sample[duplicated(tmp_sample)] <-
                2*tmp_sample[duplicated(tmp_sample)]
            }
            if (length(tmp_sample) != 1) {
              trueFc[t, tmpTraj] <- sort(tmp_sample, decreasing = T)
            } else {
              trueFc[t, tmpTraj] <- tmp_sample
            }

          }
        }
      }
    }
  }
  return(trueFc)
}

