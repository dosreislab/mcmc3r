# mcmc3r
An R package with utility functions to work with MCMCtree

You can install the package by typing in an R prompt:
```
devtools::install_github("dosreislab/mcmc3r")
```

## Plotting calibration densities
The packge can be used to plot the L and B calibration densities used by MCMCtree. Note you can use the `sn` package (available from CRAN) to plot the Skew-normal and Skew-t calibration densities. For example, suppose you have a minimum bound calibration `L(10, .1, 1, .025)`, you can plot this with:

```
curve(dL(x, 10, .1, 1, .025), n=5e2, from=0, to=100)
```

and suppose you have the joint calibration `B(5, 10, .01, .05)`, this can be plotted with:

```
curve(dB(x, 5, 10, .01, .05), n=5e2, from=0, to=20)
```

## Marginal likelihood calculation
The package can be used to prepare control files for marginal likelihood calculation with MCMCtree to, for example, select the relaxed-clock model. Marginal likelihood calculation can be carried out by the thermodynamic integration or stepping-stones methods. A tutorial is available [here](https://dosreislab.github.io/2017/10/24/marginal-likelihood-mcmc3r.html).

## Working with quantitative morphological characters
The package is also useful for working with quantitative morphological characters. It can:  

   * Generate a morphological alignment file in phylip format in which the population noise and the correlation among characters have been taken into account. This morphological alignment is suitable for Bayesian inference of divergence times with `MCMCtree`.
   * Simulate a morphological alignment, with population noise and correlations, on a phylogeny.
   * Simulate a population sample of noisy and correlated quantitative characters.
   * Generate suitable tree and control files for use in `MCMCtree`.

## References
If using the package for marginal likelihood calculation please cite:

* dos Reis et al. (2018) [Using phylogenomic data to explore the effects of relaxed clocks and calibration strategies on divergence time estimation: Primates as a test case.](https://doi.org/10.1093/sysbio/syy001) Syst. Biol., 67: 594–615.

The block bootstrap used to estimate the error of marginal likelihood estimates is described in:

* Álvarez-Carretero et al. (2022) [A species-level timeline of mammal evolution integrating phylogenomic data](https://doi.org/10.1038/s41586-021-04341-1) Nature, 602: 263–267.

The continuous morphological models are described in:

* Álvarez-Carretero et al. (2019) [Bayesian estimation of species divergence times using correlated quantitative characters.](https://doi.org/10.1093/sysbio/syz015) Syst. Biol., 68: 967–986.  
