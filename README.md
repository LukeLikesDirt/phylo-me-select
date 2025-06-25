# phyloMEM: Reducing Phylogenetic Autocorrelation via Moran's Eigenvector Maps

## Overview

`phyloMEM()` provides a phylogenetic adaptation of the MIR (minimisation of Moran’s I in the mesiduals) method implemented in [`adespatial::mem.select()`](https://github.com/adeverse/adespatial/). While `mem.select()` was designed to reduce spatial autocorrelation in ecological models, `phyloMEM()` uses a phylogenetic proximity matrix—such as those derived from patristic or Abouheif distances—to achieve the same goal in a phylogenetic context.

The method selects a subset of positive Moran’s eigenvectors (MEMs) that, when added to a model, reduce residual phylogenetic autocorrelation. The approach is especially useful when:

* You aim to test functional or ecological predictors independently of phylogenetic structure
* Your model residuals show significant phylogenetic signal, violating independence assumptions
* You are interested in removing rather than modelling the phylogenetic signal (in contrast to approaches like PGLS or `phylolm()`)

## Installation

You can install the development version of `phyloMEM` directly from GitHub using the `remotes` package:

    # Install remotes
    install.packages("remotes")

    # Install phyloMEM from GitHub
    if (!requireNamespace("remotes", quietly = TRUE)) {
      install.packages("remotes")
    }
    remotes::install_github("LukeLikesDirt/phyloMEM")

Alternatively, you can clone the repository and install manually:

    git clone https://github.com/LukeLikesDirt/phyloMEM.git
    R CMD INSTALL phyloMEM

## Key Features

* Uses [`adephylo::orthobasis.phylo()`](https://github.com/adeverse/adephylo/) to compute MEMs from a phylogenetic tree
* Uses Abouheif’s Cmean and permutation-based testing via [`adephylo::abouheif.moran()`](https://github.com/adeverse/adephylo/)
* Mirrors `mem.select()` from adespatial but is tailored for phylogenetic trees rather than spatial coordinates

## Background

This approach is based on and inspired by:

* Münkemüller, T., Lavergne, S., Bzeznik, B., Dray, S., Jombart, T., Schiffers, K., & Thuiller, W. (2012). How to measure and test phylogenetic signal. Methods in Ecology and Evolution, 3(4), 743-756
* Bauman, D., Drouet, T., Dray, S., & Vleminckx, J. (2018). Disentangling good from bad practices in the selection of spatial or phylogenetic eigenvectors. Ecography, 41(10), 1638-1649

## Example

Load required packages

    library(ape)
    library(adephylo)

Simulate a phylogenetic tree
  
    set.seed(1986)
    tree <- rtree(50)

Simulate a trait with phylogenetic signal
  
    trait <- rTraitCont(tree, model = "BM", sigma = 1)
    names(trait) <- tree$tip.label

Simulate a weakly phylogenetically structured predictor

    predictor_phylo <- rTraitCont(tree, model = "BM", sigma = 0.3)
    predictor_noise <- rnorm(length(tree$tip.label), sd = 0.7)
    predictor <- predictor_phylo + predictor_noise
    names(predictor) <- tree$tip.label

Check the names and tips match

    identical(names(trait), tree$tip.label)

Fit a base model with residual phylogenetic autocorrelation

    base_model <- lm(trait ~ predictor)
    residuals_base <- residuals(base_model)
    names(residuals_base) <- tree$tip.label

Check residual phylogenetic signal

    prox <- adephylo::proxTips(tree, method = "Abouheif")
      
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
  
    moran_resid <- phylo.moran.test(residuals_base, prox)
  
    cat("Base model residual phylogenetic signal: Abouheif’s C =", moran_resid$obs,
        ", p =", moran_resid$pvalue, "\n")

Compute phylogenetic MEMs

    ME_all <- orthobasis.phylo(tree, method = "Abouheif")
    attr(ME_all, "values")[1:5]  # First few eigenvalues

Apply MIR selection on residuals

    res <- phylo.mem.select(residuals_base, ME_all, prox, verbose = TRUE)

Add selected MEMs to base model

    if (!is.null(res$ME.select)) {
      MEM_selected <- res$ME.select
      model_df <- data.frame(
        trait = trait,
        predictor = predictor,
        MEM_selected
      )
      final_model <- lm(trait ~ predictor + ., data = model_df)
  
      # Test residuals of final model
      resid_final <- residuals(final_model)
      names(resid_final) <- tree$tip.label
      moran_final <- phylo.moran.test(resid_final, prox)
      cat("Final model residual phylogenetic signal: Abouheif’s C =", moran_final$obs,
          ", p =", moran_final$pvalue, "\n")
    }
