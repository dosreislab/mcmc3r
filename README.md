# mcmc3r
Utility functions to work with MCMCtree

## Marginal likelihood calculation
The package is currently useful for preparing control files for marginal likelihood calculation with MCMCtree to select the relaxed-clock model. A tutorial is available [here](https://dosreislab.github.io/2017/10/24/marginal-likelihood-mcmc3r.html).

If using the package please consider citing:

* dos Reis et al. (2018) [Using phylogenomic data to explore the effects of relaxed clocks and calibration strategies on divergence time estimation: Primates as a test case.](https://doi.org/10.1093/sysbio/syy001) Syst. Biol., 67: 594--615.

## Working with quantitative morphological characters
This package can also be used to work with quantitative morphological characters. The main functions are the following:  

   * Generate a morphological alignment file in `MCMCtree` format in which the population noise and the correlation among characters have been taken into account   
   * Simulate a phylogeny with noisy and correlated quantitative characters  
   * Simulate a population sample with noisy and correlated quantitative characters   
   * Generate a tree and a control file ready to use in `MCMCtree`

Once generated the file with the morphological alignment, the dating software `MCMCtree` can be used to estimate the divergence times
