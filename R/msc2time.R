#' Time-calibrate a multi-species phylogeny
#'
#' @description Calibrate an MCMC sample from a multi-species coalescent
#' phylogeny to absolute time using a fossil calibration.
#'
#' @param mcmc A data frame containing the MCMC output of a BPP A00 analysis
#' @param node.name A character vector of length one with the name of the node to which the calibration will be applied
#' @param bounds A numeric vector of length two with age bounds for the calibration
#'
#' @details
#' The function will obtain a sample of times from a uniform distribution
#' delimited by \code{bounds}. The sampled times are then used to
#' calculate the molecular rate, and then re-scale the relative
#' times (tau's) for the other nodes in \code{mcmc}. See Angelis and
#' dos Reis (2015) for details. For the BPP A00 analysis see Yang
#' (2015).
#'
#' @return
#' A data frame with a posterior sample of the calibrated times and
#' molecular rate.
#'
#' @references
#' K. Angelis and M. dos Reis (2015) \emph{The impact of ancestral
#' population size and incomplete lineage sorting on Bayesian
#' estimation of species divergence times.} Curr. Zool., 61: 874--885.
#'
#' Z. Yang (2015) \emph{The BPP program for species tree estimation
#' and species delimitation.} Curr. Zool., 61: 854--865.
#'
#' @examples
#' data(hominids)
#'
#' # Calibrate the hominid phylogeny with a fossil calibration of
#' # between 6.5 to 10 Ma for the human-chimp divergence.
#' calmsc <- msc2time.t(mcmc=hominids$mcmc, node="7humanchimp", bounds=c(6.5, 10))
#'
#' # posterior age of human-chimp (this is the same as the calibration)
#' plot(density(calmsc$t_7humanchimp, adj=.1), xlab="Time (Ma)", main="Human-chimp age")
#' rug(calmsc$t_7humanchimp)
#'
#' # posterior age of human-goriall (root of the tree)
#' plot(density(calmsc$t_5humanchimpgorillaorang, adj=.1), xlab="Time (Ma)", main="Root age")
#' rug(calmsc$t_5humanchimpgorillaorang)
#'
#' \dontrun{
#' # if you have the coda package installed, you can get the
#' # HPD intervals
#' coda::HPDinterval(coda::as.mcmc(calmsc))
#' }
#'
#' @author Mario dos Reis
#'
#' @export
msc2time.t <- function(mcmc, node.name, bounds) {
  time.name <- paste("tau_", node.name, sep="")
  i <- match(time.name, names(mcmc))
  n <- length(mcmc[,i])  # length of the MCMC sample
  ti <- grep("tau_", names(mcmc)) # find the unscaled times (tau's) in data frame
  tr <- runif(n, min=bounds[1], max=bounds[2]) # obtain randomly sampled times from calibration
  rate <- mcmc[,i] / tr  # calculate substitution rate
  tmcmc <- mcmc[,ti] / rate # calibrate times
  # rename the calibrated times
  new.names <- sub("^tau", "t", names(tmcmc))
  names(tmcmc) <- new.names

  return (cbind(tmcmc, rate)) # return a dataframe with the calibrated times
}

#' Time-calibrate a multi-species coalescent phylogeny.
#'
#' @description Calibrate an MCMC sample from a multi-species coalescent
#' phylogeny to absolute time using a rate/generation-time calibration.
#'
#' @param mcmc A data frame containing the MCMC output of a BPP A00 analysis
#' @param u.mean Numeric vector of length one with the mean for the per-generation molecular rate calibration
#' @param u.sd Numeric vector of length one with the SD for for the per-generation molecular rate calibration
#' @param g.mean Numeric vector of length one with the mean for the generation time calibration
#' @param g.sd Numeric vector of length one with the SD for the generation time calibration
#'
#' @details
#' \code{u.mean} and \code{u.sd}, and \code{g.mean} and \code{g.sd},
#' are used to construct gamma density calibrations for the
#' per-generation molecular rate and generation time respectively.
#' The gamma density with mean \eqn{m} and s.d. \eqn{s} has shape
#' \eqn{a = (m / s)^2} and scale \eqn{s = m / s^2}.
#'
#' The gamma densities are used to obtain random samples of the
#' per-generation rate and generation time. From these the molecular
#' rate per absolute time unit is calculated, and then used to
#' convert the relative times (tau's) to absolute divergence times.
#' The relative population sizes (theta's) are converted to effective
#' population sizes in number of individuals. See Yoder et al. (2016)
#' for an example with mouse lemurs. The BPP A00 analysis is
#' described in Yang (2015).
#'
#' @references
#' A. D. Yoder, C. R. Campbell, M. B. Blanco, M. dos Reis, J. U.
#' Ganzhorn, S. M. Goodman, K. E. Hunnicutt, P. A. Larsen, P. M.
#' Kappeler, R. M. Rasoloarison, J. M. Ralison, D. L. Swofford, and
#' D. W. Weisrock. (2016) \emph{Geogenetic patterns in mouse lemurs (genus
#' Microcebus) reveal the ghosts of Madagascar's forests past.}
#' Proc. Nat. Acad. Sci. USA., 113: 8049--8056.
#'
#' Z. Yang (2015) \emph{The BPP program for species tree estimation
#' and species delimitation.} Curr. Zool., 61: 854--865.
#'
#' @return
#' A data frame with a posterior sample of population sizes,
#' calibrated times, per-generation molecular rate, generation
#' time and molecular rate per time unit.
#'
#' @examples
#' data(microcebus)
#'
#' # Calibrate the Microcebus phylogeny to absoluate divergence times
#' calmsc <- msc2time.r(mcmc=microcebus$mcmc, u.mean=8.7e-9, u.sd=1.65e-9, g.mean=3.75, g.sd=0.375)
#'
#' # posterior age of the phylogeny's root (in thousans of years)
#' plot(density(calmsc$t_7OLMXRB / 1e3, adj=.1), xlab="Time (Ka)", main="Root age")
#' rug(calmsc$t_7OLMXRB / 1e3)
#'
#' # Posterior of the ancestral effective population at the root (in
#' # thousands of individuals)
#' plot(density(calmsc$Ne_7OLMXRB / 1e3, adj=.1), xlab="Ne (x 10^3 individuals)", main = "Ne at root")
#' rug(calmsc$Ne_7OLMXRB / 1e3)
#'
#' \dontrun{
#' # if you have the coda package installed, you can get the
#' # HPD intervals
#' coda::HPDinterval(coda::as.mcmc(calmsc))
#' }
#'
#' @author Mario dos Reis
#'
#' @export
msc2time.r <- function(mcmc, u.mean, u.sd, g.mean, g.sd) {
  n <- length(mcmc[,1])  # length of the MCMC sample
  ti <- grep("tau_", names(mcmc)) # find the unscaled times (tau's) in data frame
  ni <- grep("theta_", names(mcmc)) # find mutation-rate scaled population sizes in data frame

  # obtain alpha and beta re-parameterizations for the gamma
  u.a <- u.mean^2 / u.sd^2
  u.b <- u.mean / u.sd^2
  g.a <- g.mean^2 / g.sd^2
  g.b <- g.mean / g.sd^2

  # obtain random samples of u and g from gamma distribution
  u <- rgamma(n, u.a, u.b)
  g <- rgamma(n, g.a, g.b)
  rate <- u / g

  # obtain calibrated times and population sizes
  tmcmc <- mcmc
  tmcmc[,ti] <- mcmc[,ti] / rate
  tmcmc[,ni] <- mcmc[,ni] / (4 * u)

  # rename params
  new.names <- sub("^tau", "t", names(tmcmc))
  new.names <- sub("^theta", "Ne", new.names)
  names(tmcmc) <- new.names

  return(cbind(tmcmc, u, g, rate))
}

# TODO: Add a better description
#' A BPP A00 MCMC sample for a mouse lemur phylogeny
"microcebus"

# TODO: Add a better description
#' A BPP A00 MCMC sample for an hominid phylogeny
"hominids"
