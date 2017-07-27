# S4 Classes defined for thr SCsim package

# SCsimSet ====================================================================
# SCsimSet: Single Cell simulated dataset
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' The "SCsimSet" class
#'
#' S4 class used to store simulation parameters and resulting count matrices.
#'
#'@section Slots:
#'  \describe{
#'    \item{\code{sampleInfo}:}{Named list which contains sample information:
#'    total number of genes (nGenes), total number of cells (nCells), total
#'    number of cell population (nCells) and relative size / proportion of cell
#'    populations (pPop).}
#'
#'
#'    \item{\code{baseMeans}:}{Object of class \code{"DistrSet"}, which contains
#'    distribution parameters of gene basal expression.}
#'    \item{\code{cellBiais}:}{Object of class \code{"DistrSet"}, which contains
#'    distribution parameters of cell-to-cell biais (capture efficiency, batch,
#'    library size...)}
#'    \item{\code{FcDeg}:}{Object of class \code{"FCSet"}, which contains DEG
#'    information and FC estimates.}
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
#'
#'}
#'
#' @name SCsimSet
#' @rdname SCsimSet
#' @aliases SCsimSet-class
#' @references  Thanks to the scater and scran packages
#' (github.com/) for their Single Cell Expression Set class.
#' @exportClass SCsimSet

setClass("SCsimSet",
         slots = c(sampleInfo = "list",
                   baseMeans = "DistrSet",
                   cellBiais = "DistrSet",
                   FcDeg = "FCSet",
                   countDispersion = "numeric",
                   dropoutPct = "numeric",
                   dropout = "matrix",
                   effectiveMeans = "matrix",
                   baseCounts = "matrix",
                   effectiveCounts = "matrix")
)
