# mcmc3r
Utility functions to work with MCMCtree

## Marginal likelihood calculation
The package can be used to prepare control files for marginal likelihood calculation with MCMCtree to, for example, select the relaxed-clock model. Marginal likelihood calculation can be carried out by the thermodynamic integration or stepping-stones methods. A tutorial is available [here](https://dosreislab.github.io/2017/10/24/marginal-likelihood-mcmc3r.html).

## Working with quantitative morphological characters
The package is also useful for Working with quantitative morphological characters. It can:  

   * Generate a morphological alignment file in phylip format in which the population noise and the correlation among characters have been taken into account. This morphological alignment is suitable for Bayesian inference of divergence times with `MCMCtree`.
   * Simulate a morphological alignment, with population noise and correlations, on a phylogeny.
   * Simulate a population sample of noisy and correlated quantitative characters.
   * Generate suitable tree and a control files for use in `MCMCtree`.

## References
If using the package for marginal likelihood calculation please cite:

* dos Reis et al. (2018) [Using phylogenomic data to explore the effects of relaxed clocks and calibration strategies on divergence time estimation: Primates as a test case.](https://doi.org/10.1093/sysbio/syy001) Syst. Biol., 67: 594--615.

The continuous morphological models are described in:

* Álvarez-Carretero et al. (2018) [Bayesian estimation of species divergence times using correlated quantitative characters](https://doi.org/10.1101/441105) bioRxiv, 441105.
