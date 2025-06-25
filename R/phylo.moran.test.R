#' Compute Moran's I test with permutation
#'
#' Performs a permutation-based test of phylogenetic autocorrelation using Abouheif's C statistic.
#'
#' @param x A named numeric vector (e.g., trait values or residuals).
#' @param prox A phylogenetic proximity matrix (e.g., from adephylo::proxTips()).
#' @param nperm Number of permutations (default = 999).
#' @param one.tailed Logical; whether to compute a one-tailed p-value (default = TRUE).
#'
#' @return A list with:
#' \describe{
#'   \item{obs}{Observed Moran's I (Abouheif’s C)}
#'   \item{sim}{Simulated null distribution of Moran's I}
#'   \item{pvalue}{P-value for the test}
#' }
#' @export
#'
#' @examples
#' if (requireNamespace("adephylo", quietly = TRUE)) {
#'   tree <- ape::rtree(20)
#'   trait <- stats::rnorm(20)
#'   names(trait) <- tree$tip.label
#'   prox <- adephylo::proxTips(tree, method = "Abouheif")
#'   phylo.moran.test(trait, prox, nperm = 99)
#' }

phylo.moran.test <- function(x, prox, nperm = 999, one.tailed = TRUE) {
  obs_stat <- adephylo::abouheif.moran(x, prox)$obs
  sim_stats <- replicate(nperm, {
    x_perm <- sample(x)
    adephylo::abouheif.moran(x_perm, prox, nrepet = 0)$obs
  })
  p_value <- if (one.tailed) {
    (sum(sim_stats >= obs_stat) + 1) / (nperm + 1)
  } else {
    (sum(abs(sim_stats) >= abs(obs_stat)) + 1) / (nperm + 1)
  }
  return(list(obs = obs_stat, sim = sim_stats, pvalue = p_value))
}
