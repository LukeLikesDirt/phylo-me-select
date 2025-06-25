#' Phylogenetic MEM Selection via MIR to Reduce Positive Autocorrelation
#'
#' This function adapts the MIR (minimisation of Moran's I in the residuals) 
#' approach for phylogenetic data. It selects a subset of phylogenetic Moran's 
#' eigenvectors (ME) that minimise residual phylogenetic autocorrelation, 
#' focusing on positively autocorrelated eigenvectors.
#'
#' @param x A named numeric vector (e.g., trait values or model residuals).
#' @param ME_all A matrix of all Moran's eigenvectors (e.g., from orthobasis.phylo()).
#' @param prox A phylogenetic proximity matrix (e.g., from adephylo::proxTips()).
#' @param ME.all Logical; if TRUE, the complete set of ME variables is returned.
#' @param nperm Number of permutations for sequential tests (default = 999).
#' @param nperm.global Number of permutations for the global test (default = 9999).
#' @param alpha Significance threshold (default = 0.05).
#' @param verbose Logical; if TRUE, prints diagnostic messages.
#'
#' @return A list containing:
#' \describe{
#'   \item{global.test}{Global Moran's I test for residual phylogenetic autocorrelation.}
#'   \item{ME.all}{All Moran's eigenvectors (if ME.all = TRUE).}
#'   \item{ME.select}{Selected subset of positive ME variables.}
#'   \item{summary}{Summary data frame of selected variables, Moran's I, and p-values.}
#' }
#'
#' @export
#'
#' @examples
#' if (requireNamespace("ape", quietly = TRUE) && requireNamespace("adephylo", quietly = TRUE)) {
#'   set.seed(123)
#'   tree <- ape::rtree(30)
#'   trait <- ape::rTraitCont(tree, model = "BM")
#'   names(trait) <- tree$tip.label
#'
#'   # Compute MEMs and proximity matrix
#'   ME_all <- adephylo::orthobasis.phylo(tree, method = "Abouheif")
#'   prox <- adephylo::proxTips(tree, method = "Abouheif")
#'
#'   # Run MEM selection
#'   result <- phylo.mem.select(trait, ME_all, prox, verbose = TRUE)
#'
#'   # View summary
#'   print(result$summary)
#' }
