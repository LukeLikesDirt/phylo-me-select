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
#'

phylo.mem.select <- function(
    x, ME_all, prox,
    ME.all = FALSE, nperm = 999, nperm.global = 9999,
    alpha = 0.05, verbose = FALSE
) {
  if (any(is.na(x))) stop("NA entries in x")
  if (length(x) != nrow(ME_all)) stop("Length of x must match rows in ME_all")
  if (is.null(names(x))) stop("x must be a named vector")
  if (!all(names(x) %in% rownames(ME_all))) stop("Names in x must match rownames in ME_all")

  x <- x[rownames(ME_all)]
  moran_values <- attr(ME_all, "values")
  positive_idx <- which(moran_values > 0)

  if (length(positive_idx) == 0) {
    if (verbose) message("No eigenvectors with positive phylogenetic autocorrelation found")
    result <- list(global.test = list(obs = NA, pvalue = 1))
    if (ME.all) result$ME.all <- ME_all
    return(result)
  }

  ME_positive <- ME_all[, positive_idx, drop = FALSE]
  positive_moran <- moran_values[positive_idx]
  global_moran <- phylo.moran.test(x, prox, nperm = nperm.global)

  res <- if (ME.all) list(global.test = global_moran, ME.all = ME_all) else list(global.test = global_moran)
  if (global_moran$pvalue >= alpha) {
    if (verbose) message("No significant positive phylogenetic autocorrelation")
    return(res)
  }

  if (verbose) message(sprintf("Global test significant (p = %.4f). Starting MIR selection with %d ME variables.",
                               global_moran$pvalue, ncol(ME_positive)))

  p <- global_moran$pvalue
  ME.sel <- idx.min <- min.moran <- p.vector <- c()
  current_x <- x

  while (p < alpha) {
    I.vector <- sapply(1:ncol(ME_positive), function(i) {
      if (i %in% ME.sel) return(NA)
      temp_lm <- lm(current_x ~ ME_positive[, i])
      temp_resid <- residuals(temp_lm)
      names(temp_resid) <- names(current_x)
      phylo.moran.test(temp_resid, prox, nperm = 1)$obs
    })

    idx.min <- which.min(abs(I.vector))

    if (verbose) message(sprintf("Testing ME variable %d (original index: %d)",
                                 length(ME.sel) + 1, positive_idx[idx.min]))

    current_x <- residuals(lm(current_x ~ ME_positive[, idx.min]))
    names(current_x) <- names(x)

    testI <- phylo.moran.test(current_x, prox, nperm = nperm)
    p <- testI$pvalue

    ME.sel <- c(ME.sel, idx.min)
    min.moran <- c(min.moran, testI$obs)
    p.vector <- c(p.vector, p)
  }

  if (verbose) message(sprintf("Procedure stopped (p = %.4f > alpha = %.2f)", p, alpha))

  selected_ME <- ME_positive[, ME.sel, drop = FALSE]
  original_indices <- positive_idx[ME.sel]

  summary_df <- data.frame(
    variables = paste0("ME", original_indices),
    order = original_indices,
    moran_eigenvalue = positive_moran[ME.sel],
    Iresid = min.moran,
    pvalue = p.vector,
    row.names = 1:length(ME.sel)
  )

  res$ME.select <- selected_ME
  res$summary <- summary_df
  return(res)
}

