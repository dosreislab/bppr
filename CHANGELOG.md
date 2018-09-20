# Changelog
Important changes to this project will be documented in this file.

We try to follow [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and we use [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.5.0] - hehehe
### Added
- Function `mcmc.summary` to generate a summary (mean, CIs) from an MCMC sample from BPP or MCMCTree.
- Function `mcmc2densitree` to plot a densi-tree using an MCMC sample from BPP or MCMCTree.

### Fixed
- A bug in `mcmc2multiphylo` where the MCMC sample was not being thinned correctly.

## [0.4.0] - 2018-09-08
### Added
- Function `mcmc2multiphylo` to convert MCMC samples from BPP or MCMCTree to lists of trees.

**Tags:** Added, Changed, Deprecated, Removed, Fixed, Security
