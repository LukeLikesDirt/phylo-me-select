# `phyloMEM`: Reducing Phylogenetic Autocorrelation via Moran's Eigenvector Maps

## Overview

`phyloMEM` provides a phylogenetic adaptation of the MIR (minimisation of Moran’s I in the mesiduals) method implemented in `adespatial::mem.select()`. While `mem.select()` was designed to reduce spatial autocorrelation in ecological models, `phyloMEM` uses a phylogenetic proximity matrix—such as those derived from patristic or Abouheif distances—to achieve the same goal in a phylogenetic context.

The method selects a subset of positive Moran’s eigenvectors (MEMs) that, when added to a model, reduce residual phylogenetic autocorrelation. The approach is especially useful when:

    You aim to test functional or ecological predictors independently of phylogenetic structure.

    Your model residuals show significant phylogenetic signal, violating independence assumptions.

    You are interested in removing rather than modelling the phylogenetic signal (in contrast to approaches like PGLS or `phylolm()`).

## Key Features

    Uses `adephylo::orthobasis.phylo()` to compute MEMs from a phylogenetic tree.

    Uses Abouheif’s Cmean and permutation-based testing via `adephylo::abouheif.moran()`.

    Mirrors `mem.select()` from adespatial but is tailored for phylogenetic trees rather than spatial coordinates.

## Background

This approach is based on and inspired by:

    Dray, S., Legendre, P., & Peres-Neto, P. R. (2006). Spatial modelling: a comprehensive framework for principal coordinate analysis of neighbour matrices (PCNM). Ecological Modelling, 196(3–4), 483–493.

    Dray, S., Bauman, D., Blanchet, F. G., Borcard, D., Clappe, S., Guénard, G., ... & Wagner, H. H. (2022). adespatial: Multivariate Multiscale Spatial Analysis. R package version 0.3-21.

    Pavoine, S., Ollier, S., & Dufour, A. B. (2005). Is the originality of a species measurable? Ecology Letters, 8(6), 579–586.
