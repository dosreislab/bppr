# bppr: a helper R package for BPP

Currently the package is useful for

* Calibrating BPP phylogenies to geological time

* Calculation of Bayes factors with BPP for model selection

The package can calibrate a BPP A00 analysis to absolute divergence times by using either a fossil calibration on a node age, or a prior on the per-generation rate and generation time. In the latter case a posterior sample of the effective population sizes is obtained.

Bayes factor calculations are useful for species delimitation with very large datasets, in which case the rjMCMC algorithm may be inefficient. Bayes factors with bppr are calculated with the stepping stones algorithm or the Gaussian quadrature (thermodynamic integration) approach of Rannala and Yang (2018). Note that the stepping stones algorithm appears to be much more efficient than the Gaussian quadrature method.

A tutorial for the package can be found [here](https://dosreislab.github.io/2018/08/31/bppr.html).

## Installation

If you have the devtools package installed, you can install bppr by typing in R:

```R
devtools::install_github("dosreislab/bppr")
```

## Example
Calibrating the hominid phylogeny to geological time and plotting it:

```R
data(hominids)
# Calibrate the hominid phylogeny with a uniform fossil calibration of
# between 6.5 to 10 Ma for the human-chimp divergence, and plot the
# calibrated sample
calmsc <- msc2time.t(mcmc=hominids$mcmc, node="7humanchimp", calf=runif,
  min=6.5, max=10)
mcmc2densitree(hominids$tree, calmsc, "t_", thin=0.05, alpha=0.01)
  title(xlab="Divergence time (Ma)")
```

![](/figs/apes.png)

## References

If using this package, please consider citing the following whether appropriate.

The A00 analysis is described in:

* Z. Yang (2015) [The BPP program for species tree estimation and species delimitation](https://doi.org/10.1093/czoolo/61.5.854). Curr. Zool., 61: 854--865.

The general random sampling procedure to recalibrate the divergence times is given in:

* K. Angelis and M. dos Reis (2015) [The impact of ancestral population size and incomplete lineage sorting on Bayesian estimation of species divergence times](https://doi.org/10.1093/czoolo/61.5.874). Curr. Zool., 61: 874--885.

The random sampling procedure using the per-generation rate and generation time is described in:

* A. D. Yoder et al. (2016) [Geogenetic patterns in mouse lemurs (genus Microcebus) reveal the ghosts of Madagascar's forests past](https://doi.org/10.1073/pnas.1601081113). Proc. Nat. Acad. Sci. USA., 113: 8049--8056.

Gaussian quadrature for Bayes factors in BPP is described in:

* Rannala, B., and Z. Yang. 2017. [Efficient Bayesian species tree inference under the multispecies coalescent](https://doi.org/10.1093/sysbio/syw119). Systematic Biology, 66: 823-842.

Other relevant citations are given in the helpfiles of the package.
