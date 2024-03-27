# Changelog
Important changes to this project will be documented in this file.

We try to follow [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and we use [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.5.5] - 2024-03-27

### Removed

- Update `DESCRIPTION` file:
  - Remove `Morpho` and `Rdpack` from `Import`.
  
### Fixed

- Update `Morpho.R` script:
  - Check whether users have `Morpho` installed when running function
   `proc2MCMCtree`, stop if not -> error and warning solved.
  - Update `.checkData` and `.checkArr` functions so that checks with `class()`
    do not raise warnings -> warning solved.
- Update `calibrations.R` file:
  - Remove additional parenthesis in line 45 (i.e. from `pv <- pL(qv, tL=1))`
    to `pv <- pL(qv, tL=1))`) -> error solved.
- Update `dBD.R` script:
  - Define parameter `x` as `@param` -> warning solved.
- Update `mcmc2densitree`:
  - Change `cex.labels` to `cex.lab` (line 18, `@param` section) as the latter is
    used in the function -> warning solved.

### Changed

- Update `DESCRIPTION` file:
  - Upgrade `RoxygenNote` to 7.3.1.
  - Upgrade minimum R version to `(>= 3.5.0)`.
- Update `Morpho.R` script:
  - Add `\dontrun{}` for examples in `proc2MCMCtree` that require package
    `Morpho` in case users do not have this package installed.
- Update `mcmc2anc.R` script:
  - Add `\dontrun{}` for examples section to avoid warnings when using `rgl`.
  - Change lines 45-46 so that lines are not wider than 100 characters as
    before they were truncated in the PDF manual.
- Update vignette `Reproduce_Carnivora_analysis.Rmd`:
  - Add `eval=FALSE` in all code snippets. Package `Morpho` is not imported and
    running these examples raises errors.
  - Define code snippets from section 4 onward using `{r}`.
  - Add `mcm3r::` before calling a function within this package in the code
    snippets.
  - Replace some `=` with `<-` as it is encouraged in R.
  - Fix some typos.
- Update `carnivores.R` script:
  - Fix typo in surname.

## [0.5.4] - 2024-02-09
### Added
- Function `mcmc.sum` to calculate summary statistics on MCMCtree's output.

## [0.5.3] - 2023-12-07
### Added
- Functions `pL` and `qL` to calculate the probability and quantile functions
for the minim bound (truncated-Cauchy) calibrations.

## [0.5.2] - 2023-04-24
### Added
- Function `mcmc2anc` for ancestral character reconstruction from an MCMC sample.
- Directory `misc/ape4s`, and moved all ape4s files into it.

## Changed
- The analysis of carnivores vignette, `Reproduce_Carnivora_analysis.Rmd`, to 
streamline the text and add an example of ancestral reconstruction.

## [0.5.1] - 2023-03-30
### Added
- Function `mcmc2densitree` to plot densitrees using the MCMC output from either
MCMCtree or BPP.
- `microcebus` and `hominid` BPP A00 MCMC datasets.
- Carnivora data to `misc/carnivora`, suitable for reproducing MCMC sampling 
of divergence times using molecular and morphometric data with MCMCtree.

### Changed
- The analysis of carnivores vignette, `Reproduce_Carnivora_analysis.Rmd`, to 
simplify the text and to add MCMC sampling using MCMCtree. Added files 
`figtree.jpg` and `densitree.jpg` to illustrate the vignette.

## [0.4.7] - 2023-02-02
### Fixed
- An error was returned in several instances of testing input where class(x) 
was used to check matrices due to a new behaviour of classes in matrices from 
R > v4.x. Now `inherits` is used to test matrix classes instead of `class`.

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
- Function `dBL` to calculate the kernel density for the birth-death process
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
