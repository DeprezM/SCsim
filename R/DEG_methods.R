# describeDEG =================================================================

#' describeDEG -> Differentially expressed genes indexes / type / regulation in
#' the gene expression table to generate
#'
#'
#'
#'
#'
#'
#'
#'
#'


DEG_table <- describeDEG(nGenes = 1000, nCells = 100, nPop = 3,
                         pPop = c(15, 35, 50), pDEG = c(15, 35, 50),
                         pDE = 40, pDP = 20, pDM = 30, pDC = 10,
                         pUp = 50, pDown= 50)



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


    if (length(pDE) > 1 & length(pDM) > 1){
      for (i in 1:nPop){
        if (!isTRUE(all.equal(pDC + pDE[i] + pDP + pDM[i], 100))){
        stop("sum of pDC, pDE, pDP and  pDM != 100. Those proportions are
             sub-types of DEG which sum must be equal to 100.")
        }
      }
    } else if (!isTRUE(all.equal(pDC + pDE + pDP + pDM, 100))){
      stop("sum of pDC, pDE, pDP and  pDM != 100. Those proportions are
             sub-types of DEG which sum must be equal to 100.")
    } else {
      pDE <- rep(pDE, nPop)
      pDM <- rep(pDM, nPop)
    }

    # Check percentage validity
    if (!isTRUE(all.equal(pUp + pDown, 100))) {
      stop("sum of pUp and  pDown != 100.")
    } else {
      pUp <- pUp / 100
      pDown <- pDown / 100
    }

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
      print(pop)
      genesStatus$population[popIndex:(popIndex + nDEGpop[pop] - 1)] <-
        rep(pop, nDEGpop[pop])
      print(nDEGpop[pop])
      geneType <- c()
      geneReg <- c()

      if ("pDE" %in% listDEGtype) {
        if (pop == 1) {
          pDE <- pDE/100
        }
        geneType <- c(geneType, rep("DE", ceiling(nDEGpop[pop] * pDE[pop])))
        geneReg <- c(geneReg, rep("Up",
                              ceiling(ceiling(nDEGpop[pop] * pDE[pop]) * pUp)))
        geneReg <- c(geneReg, rep("Down", sum(geneType == "DE") -
                              ceiling(ceiling(nDEGpop[pop] * pDE[pop]) * pUp)))
      }

      if ("pDM" %in% listDEGtype) {
        if (pop == 1) {
          pDM <- pDM/100
        }
        geneType <- c(geneType, rep("DM", ceiling(nDEGpop[pop] * pDM[pop])))
        geneReg <- c(geneReg, rep("Up",
                               ceiling(ceiling(nDEGpop[pop] * pDM[pop]) * pUp)))
        geneReg <- c(geneReg, rep("Down", sum(geneType == "DM") -
                               ceiling(ceiling(nDEGpop[pop] * pDM[pop]) * pUp)))

      }

      if ("pDP" %in% listDEGtype) {
        if (pop == 1) {
          pDP <- pDP/100
        }
        geneType <- c(geneType, rep("DP", ceiling(nDEGpop[pop] * pDP)))
        geneReg <- c(geneReg, rep("Up",
                                  ceiling(ceiling(nDEGpop[pop] * pDP) * pUp)))
        geneReg <- c(geneReg, rep("Down", sum(geneType == "DP") -
                                    ceiling(ceiling(nDEGpop[pop] * pDP) * pUp)))
      }
      if ("pDC" %in% listDEGtype) {
        if (pop == 1) {
          pDC <- pDC/100
        }
        geneType <- c(geneType, rep("DC", ceiling(nDEGpop[pop] * pDC)))
        geneReg <- c(geneReg, rep("Up",
                                  ceiling(ceiling(nDEGpop[pop] * pDC) * pUp)))
        geneReg <- c(geneReg, rep("Down", sum(geneType == "DC") -
                                    ceiling(ceiling(nDEGpop[pop] * pDC) * pUp)))
      }
      geneType <- c(rep("DE", nDEGpop[pop]-length(geneType)), geneType)
      geneReg <- c(rep("Up", nDEGpop[pop]-length(geneType)), geneReg)

      print(length(geneType))
      print(length(geneReg))

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

