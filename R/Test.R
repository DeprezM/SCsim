# Test parameters

nGenes = 1000
nCells = 200
nPop = 3
pPop = c(30, 10, 60)
seed = 25

distribution = "gamma"
fstParam = 2
sndParam = 0.75

nbBatch = 2
cellsPerBatch = c(30, 70)
batchEffect = NULL

pDEG = c(20, 30, 10)
pDE = c(40, 10, 20)
pDP = 20
pDM = c(30, 60, 50)
pDC = 10
pUp = c(70, 50, 20)
pDown= c(30, 50, 80)

cellMixedDP = "pseudo"
mixDP = 25
cellMixedDM = "pseudo"
mixDM = 25
popMixDP = NULL
trajectory = list(c(1,2,3))
doublet = 2

distrUpFc = "medium"
distrDownFc = "medium"

dropoutPct = 50


# t <- newSCsimSet(nGenes = 1000,
#                  nCells = 200,
#                  nPop = 3,
#                  pPop = c(30, 10, 60),
#                  seed = 25,
#
#                  distribution = "gamma",
#                  fstParam = 2,
#                  sndParam = 0.75,
#
#                  nbBatch = 2,
#                  cellsPerBatch = c(30, 70),
#                  batchEffect = NULL,
#
#                  pDEG = c(20, 30, 10),
#                  pDE = c(40, 10, 20),
#                  pDP = 20,
#                  pDM = c(30, 60, 50),
#                  pDC = 10,
#                  pUp = c(70, 50, 20),
#                  pDown= c(30, 50, 80),
#
#                  cellMixedDP = "pseudo",
#                  mixDP = 25,
#                  cellMixedDM = "pseudo",
#                  mixDM = 25,
#                  popMixDP = NULL,
#                  trajectory = list(c(1,2,3)),
#                  doublet = 2,
#
#                  distrUpFc = "medium",
#                  distrDownFc = "medium",
#
#                  dropoutPct = 50
# )

