# Changelog
Important changes to this project will be documented in this file.

We try to follow [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and we use [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.6.3] - 2023-03-28
### Fixed
- Function `mcmc2densitree` so that tip labels are correctly positioned in the
plot when extinct taxa are present (set using tip.ages).

## [0.6.2] - 2022-05-25
### Added
- Argument `tip.ages` to `mcmc2densitree`. This allows plotting trees with 
sequences that have different sampling times (such as viruses).

## [0.6.1] - 2020-11-29
### Added
- An additional example to the `msc2time` functions documentation using the
shifted log-normal distribution as a calibration density.

## [0.6.0] - 2020-02-25
### Changed
- Function `bayes.factors` so that it calculates parameteric bootstrap of
posterior probabilities. Updated documentation.

## [0.5.3] - 2019-11-12
### Fixed
- A bug where the Gaussian quadrature rules (`glqrules`) could not be found.
This meant method "gauss-quad" could not be used with the `make.beta` function.

## [0.5.2] - 2019-01-17
### Changed
- Small changes to examples in `msc2time` and `mcmc.summary`

## [0.5.1] - 2018-09-20
### Added
- An additional example to the `mcmc2densitree` helpfile

## [0.5.0] - 2018-09-20
### Added
- Function `mcmc.summary` to generate a summary (mean, CIs) from an MCMC sample
from BPP or MCMCTree.
- Function `mcmc2densitree` to plot a densi-tree using an MCMC sample from BPP
or MCMCTree.

### Fixed
- A bug in `mcmc2multiphylo` where the MCMC sample was not being thinned
correctly.

## [0.4.0] - 2018-09-08
### Added
- Function `mcmc2multiphylo` to convert MCMC samples from BPP or MCMCTree to
lists of trees.

**Tags:** Added, Changed, Deprecated, Removed, Fixed, Security
