% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc.summary.R
\name{mcmc.summary}
\alias{mcmc.summary}
\title{MCMC summary}
\usage{
mcmc.summary(mcmc, prob = 0.95)
}
\arguments{
\item{mcmc}{Data frame with the MCMC output of BPP or MCMCTree}

\item{prob}{Numeric, probability for credibility interval calculation}
}
\value{
A list with elements \code{means}, \code{eq.ci}, and \code{hpd.ci}
  containing the posterior (or prior) means, equal tail CI and HPD CI.
}
\description{
Create an MCMC summary from a BPP or MCMCTree analysis
}
\details{
\code{mcmc} should contain the output (say from file
  \code{mcmc.txt}) generated by a BPP A00 or MCMCtree analysis. The function
  will calculate the posterior (prior) means, the equal-tail credibility
  interval (CI), and the highest posterior (prior) density (HPD) CI. \code{prob} is
  used to calculate the CIs. For example, if \code{prob} = 95%, then you get the
  95% CIs. This function requires the CODA pacakge to calculate the HPD CIs.
}
\examples{
\dontrun{
mcmc.summary(hominids$mcmc[,-1])
}

}
\author{
Mario dos Reis
}
