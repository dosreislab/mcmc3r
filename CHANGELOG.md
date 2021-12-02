# Changelog
Important changes to this project will be documented in this file.

We try to follow [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and we use [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.4.6] - 2021-12-02
### Fixed
- A bug in the `dBD` example (missing comment character). 

## [0.4.5] - 2021-07-15
### Changed
- Function `gauss.quad` so it now accepts a `se` argument (as for 
`stepping.stones`). This solves a bug when running `gauss.quad.boot`.

## [0.4.4] - 2021-07-14
### Added
- Function `gauss.quad.boot` to calculate marginal likelihood on bootstrap
replicates using the Gaussian quadrature method. Documentation for
`gauss.quad.boot` and `stepping.stones.boot` are merged.

## [0.4.3] - 2020-03-24
### Added
- Function 'dBL' to calculate the kernel density for the birth-death process
with species sampling.

## [0.4.2] - 2020-03-13
### Added
- Functions `dL`, `dB` and `dU` to calcute the fossil calibration densites used
in MCMCtree (useful for plotting the densities).

## [0.4.1] - 2020-02-24
### Fixed
- An error in `block.boot` that caused generation of empty likelihood block 
samples.

## [0.4.0] - 2020-12-31
### Added 
- Functions `block.boot` and `stepping.stones.boot` to perform block (stationary)
bootstrap resampling to estimate the standard error for the log-marginal 
likelihood under the stepping stones method. Suitable for short MCMC runs for
which the delta approximation to estimate the standard error may not work well.

## [0.3.2] - 2020-12-30
### Fixed
- Functions `stepping.stones` and `gauss.quad` so that they work correctly 
with incomplete mcmctree MCMC output files.

## [0.3.1] - 2020-02-24
### Changed
- Function `bayes.factors` so that bootstrap CI's are printed with models as
rows and 95% CI as columns.

## [0.3.0] - 2018-10-14
### Added
A collection of functions for working with continuous morphological data. These
are the functions used in Álvarez-Carretero et al. (2018, bioRxiv, 441105).

## [0.2.0] - 2018-10-08
### Changed
- Function `bayes.factors` so that it now calculates bootstrap CIs for posterior
probabilities. Prior probabilities can now be used.

**Tags:** Added, Changed, Deprecated, Removed, Fixed, Security
