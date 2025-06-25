# phyloMEM: Reducing Phylogenetic Autocorrelation via Moran's Eigenvector Maps

## Overview

`phyloMEM()` provides a phylogenetic adaptation of the MIR (minimisation of Moran’s I in the mesiduals) method implemented in [`adespatial::mem.select()`](https://github.com/adeverse/adespatial/). While `mem.select()` was designed to reduce spatial autocorrelation in ecological models, `phyloMEM()` uses a phylogenetic proximity matrix—such as those derived from patristic or Abouheif distances—to achieve the same goal in a phylogenetic context.

The method selects a subset of positive Moran’s eigenvectors (MEMs) that, when added to a model, reduce residual phylogenetic autocorrelation. The approach is especially useful when:

* You aim to test functional or ecological predictors independently of phylogenetic structure
* Your model residuals show significant phylogenetic signal, violating independence assumptions
* You are interested in removing rather than modelling the phylogenetic signal (in contrast to approaches like PGLS or `phylolm()`)

## Key Features

* Uses [`adephylo::orthobasis.phylo()`](https://github.com/adeverse/adephylo/) to compute MEMs from a phylogenetic tree
* Uses Abouheif’s Cmean and permutation-based testing via [`adephylo::abouheif.moran()`](https://github.com/adeverse/adephylo/)
* Mirrors `mem.select()` from adespatial but is tailored for phylogenetic trees rather than spatial coordinates

## Background

This approach is based on and inspired by:

* Münkemüller, T., Lavergne, S., Bzeznik, B., Dray, S., Jombart, T., Schiffers, K., & Thuiller, W. (2012). How to measure and test phylogenetic signal. Methods in Ecology and Evolution, 3(4), 743-756
* Bauman, D., Drouet, T., Dray, S., & Vleminckx, J. (2018). Disentangling good from bad practices in the selection of spatial or phylogenetic eigenvectors. Ecography, 41(10), 1638-1649

## Example

