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
#' The MCMC sample in \code{mcmc} must have been obtained under a BPP A00
#' analysis.
#'
#' The function will obtain a sample of times from a uniform distribution
#' delimited by \code{bounds}. The sampled times are then used to
#' calculate the molecular rate, and then re-scale the relative
#' times (tau's) for the other nodes in \code{mcmc}.
#'
#' @return
#' A data frame with the calibrated times and molecular rate.
#'
#' @examples
#' data(hominids)
#'
#' # Calibrate the hominid phylogeny with a fossil calibration of
#' # between 6.5 to 10 Ma for the human-chimp divergence.
#' calmsc <- msc2time.t(mcmc=hominids$mcmc, node="7humanchimp", bounds=c(6.5, 10))
#'
#' # posterior age of human-chimp (this is the calibration)
#' plot(density(calmsc$t_7humanchimp, adj=.1), xlab="Time (Ma)")
#' rug(calmsc$t_7humanchimp)
#'
#' # posterior age of human-goriall (root of the tree)
#' plot(density(calmsc$t_5humanchimpgorillaorang, adj=.1), xlab="Time (Ma)")
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

#' Calibrate a BPP A00 analysis using a prior on the rate
#' @author Mario dos Reis
#'
#' @export
msc2time.r <- function(mcmc, u.mean, u.sd, g.mean, g.sd) {
  n <- length(mcmc[,i])  # length of the MCMC sample
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
  tmcmc <- mcmc[,ti] / rate
  tmcmc <- mcmc[,ni] / (4 * u)

  return(cbind(tmcmc, u, g, rate))
}

# TODO: Add a better description
#' A BPP A00 MCMC sample for a mouse lemur phylogeny
"microcebus"

# TODO: Add a better description
#' A BPP A00 MCMC sample for an hominid phylogeny
"hominids"
